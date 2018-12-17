package workshop.client.manager;

import workshop.client.common.SheduleTableRenderer;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.WindowConstants;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import workshop.client.common.IGUI;
import workshop.common.Constants;

public class ManagerGUI implements IGUI {

    private static final int WIDTH = 500;
    private static final int HEIGHT = 500;

    private static final String TITLE = "Клиент менеджера";
    private static final String EDIT_ORDER = "Изменить статус";
    private static final String REFRESH = "Обновить расписание";
    private static final String DELETE = "Удалить заказ";
    private static final String CHANGE_TIME = "Изменить время";
    
    public final JFrame mFrame;
    private final JButton mRefreshButton;
    private final JButton mEditButton;
    private final JTable mSheduleTable;
    private final JScrollPane mScrollPane;
    private final JLabel mInfoLabel;
    private final JButton mDeleteButton;
    private final JButton mChangeTimeButton;
    
    public final JFrame mStatusEditFrame;
    private final JList mStatusEditList;
    private final JButton mStatusOkButton;
    private final JButton mStatusCancelButton;
    private final JTextArea mStatusDescriptionText;

    private final ClickListener mClickListener;

    private final IMController mController;
    
    private boolean isTimeChanging = false;
    private String oldTime = "";
    
    public ManagerGUI(IMController controller) {
        mController = controller;

        mFrame = new JFrame(TITLE);
        mFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        mFrame.setSize(WIDTH, HEIGHT);
        mFrame.setLayout(null);

        mClickListener = new ClickListener();

        mEditButton = new JButton(EDIT_ORDER);
        mEditButton.setBounds(50, 200, 200, 30);
        mEditButton.addActionListener(mClickListener);

        mRefreshButton = new JButton(REFRESH);
        mRefreshButton.setBounds(250, 200, 200, 30);
        mRefreshButton.addActionListener(mClickListener);

        mDeleteButton = new JButton(DELETE);
        mDeleteButton.setBounds(50, 230, 200, 30);
        mDeleteButton.addActionListener(mClickListener);

        mChangeTimeButton = new JButton(CHANGE_TIME);
        mChangeTimeButton.setBounds(250, 230, 200, 30);
        mChangeTimeButton.addActionListener(mClickListener);

        mSheduleTable = new JTable();
        mSheduleTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        mSheduleTable.setBounds(50, 50, 400, 150);

        mInfoLabel = new JLabel("");
        mInfoLabel.setBounds(50, 300, 400, 100);

        mSheduleTable.addMouseListener(new java.awt.event.MouseAdapter() {
            @Override
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                int row = mSheduleTable.rowAtPoint(evt.getPoint());
                int col = mSheduleTable.columnAtPoint(evt.getPoint());
                if (row >= 0 && col >= 0) {
                    if(!isTimeChanging) {
                    String[] info = mController.getInfo(mSheduleTable.getModel().getValueAt(row, col) +
                            "-" + mSheduleTable.getColumnName(col)).split(";");
                    if(info.length > 1) {
                        info[0] = "<html>Имя клиента: " + info[0];
                        info[1] = "<br>Описание заказа: " + info[1];
                        info[2] = "<br>Статус заказа: " + info[2];
                        info[3] = "<br>Описание статуса: " + info[3];
                        info[4] = "<br>Телефон: " + info[4] + "</html>";
                    }
                    String infoResult = Arrays.toString(info);
                    infoResult = infoResult.substring(1, infoResult.length() - 1);
                    mInfoLabel.setText(infoResult);
                    mFrame.repaint();
                    } else {
                        if (col > 9) {
                            String newTime = mSheduleTable.getModel().getValueAt(row, col)
                                    + "-" + mSheduleTable.getColumnName(col);
                            String result = mController.changeTime(oldTime, newTime);
                            isTimeChanging = false;
                            mChangeTimeButton.setText("Изменить время");
                            fillTable();
                            JOptionPane.showMessageDialog(mFrame, result);
                        } else {
                            JOptionPane.showMessageDialog(mFrame, 
                                    "Вы должны выбрать время, начиная с сегодняшнего!");
                        }
                    }
                }
            }
        });

        mScrollPane = new JScrollPane(mSheduleTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        mScrollPane.setBounds(50, 50, 400, 150);

        mStatusEditFrame = new JFrame(TITLE);
        mStatusEditFrame.setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
        mStatusEditFrame.setSize(250, 250);
        mStatusEditFrame.setLayout(null);
        
        mStatusEditList = new JList(new String[]{Constants.STATUS_WAIT_FOR_CAR,
                                                 Constants.STATUS_IN_WORK,
                                                 Constants.STATUS_NEED_FOR_DETAILS,
                                                 Constants.STATUS_REPAIRS_COMPLETED});
        mStatusEditList.setBounds(25, 25, 200, 75);
        
        mStatusOkButton = new JButton("ОК");
        mStatusOkButton.setBounds(25, 175, 100, 25);
        mStatusOkButton.addActionListener(mClickListener);
        
        mStatusCancelButton = new JButton("Отменить");
        mStatusCancelButton.setBounds(125, 175, 100, 25);
        mStatusCancelButton.addActionListener(mClickListener);
        
        mStatusDescriptionText = new JTextArea("");
        mStatusDescriptionText.setBounds(25, 110, 200, 50);
    }

    @Override
    public void show() {
        mFrame.setVisible(true);
        showStartScreen();
    }

    public void showStartScreen() {
        clear();
        mFrame.add(mRefreshButton);
        mFrame.add(mEditButton);
        mFrame.add(mDeleteButton);
        mFrame.add(mChangeTimeButton);
        mFrame.add(mScrollPane);
        mFrame.add(mInfoLabel);
        fillTable();
        mFrame.repaint();
    }
    
    public void showEditStatusScreen() {
        mStatusEditFrame.setVisible(true);
        mStatusEditFrame.add(mStatusEditList);
        mStatusEditFrame.add(mStatusOkButton);
        mStatusEditFrame.add(mStatusCancelButton);
        mStatusEditFrame.add(mStatusDescriptionText);
    }
    
    public void clear() {
        mFrame.remove(mSheduleTable);
        mFrame.remove(mScrollPane);
        mFrame.repaint();
    }

    public void clearTable() {
        DefaultTableModel tableModel = (DefaultTableModel) mSheduleTable.getModel();
        tableModel.setRowCount(0);
        tableModel.setColumnCount(0);
    }

    public void fillTable() {
        clearTable();
        DateFormat dateFormat = new SimpleDateFormat("dd.MM.yyyy");
        Calendar calendar = Calendar.getInstance();
        calendar.add(Calendar.DATE, -10);
        String currentDate = dateFormat.format(calendar.getTime());
        Object[] columns = new Object[40];
        for (int i = 0; i < 40; i++) {
            columns[i] = (String) currentDate;
            calendar.add(Calendar.DATE, 1);
            currentDate = dateFormat.format(calendar.getTime());
        }
        TableModel tableModel = new DefaultTableModel(columns, 0);
        mSheduleTable.setModel(tableModel);
        DefaultTableModel model = (DefaultTableModel) mSheduleTable.getModel();
        int hour = 10;
        for (int i = 0; i < 12; i++) {
            Object[] row = new Object[40];
            for (int j = 0; j < 40; j++) {
                row[j] = hour + ":00";
            }
            hour++;
            model.addRow(row);
        }
        mSheduleTable.setDefaultEditor(Object.class, null);

        String[] shedule = mController.getShedule();
        if (shedule == null) {
            JOptionPane.showMessageDialog(mFrame,
                    "Не удалось получить расписание.");
        } else {
            int size = shedule.length;
            Point[] filledCells = new Point[size];
            for (int i = 0; i < size; i++) {
                String[] splited = shedule[i].split("-");
                if (splited.length == 2) {
                    String time = splited[0];
                    String day = splited[1];
                    int colIndex = getColumnIndex(mSheduleTable, day);
                    if (colIndex != -1) {
                        int rowIndex = getRowIndex(time);
                        if (rowIndex != -1) {
                            filledCells[i] = new Point(rowIndex, colIndex);
                        }
                    }
                }
            }
            setTableRenderer(size, filledCells);
        }
    }
    
    private void setTableRenderer(int size, Point[] filledCells) {
        mSheduleTable.setDefaultRenderer(mSheduleTable.getColumnClass(0), new SheduleTableRenderer(size, filledCells));
    }

    private int getColumnIndex(JTable table, String header) {
        for (int i = 0; i < table.getColumnCount(); i++) {
            if (table.getColumnName(i).equals(header)) {
                return i;
            }
        }
        return -1;
    }

    private int getRowIndex(String header) {
        switch (header) {
            case "10:00": {
                return 0;
            }
            case "11:00": {
                return 1;
            }
            case "12:00": {
                return 2;
            }
            case "13:00": {
                return 3;
            }
            case "14:00": {
                return 4;
            }
            case "15:00": {
                return 5;
            }
            case "16:00": {
                return 6;
            }
            case "17:00": {
                return 7;
            }
            case "18:00": {
                return 8;
            }
            case "19:00": {
                return 9;
            }
            case "20:00": {
                return 10;
            }
            case "21:00": {
                return 11;
            }
            default: {
                return -1;
            }
        }
    }

    private class ClickListener implements ActionListener {

        @Override
        public void actionPerformed(ActionEvent e) {
            if (e.getSource() == mEditButton) {
                if(mSheduleTable.getSelectedColumn() != -1 && mSheduleTable.getSelectedRow() != -1)
                    showEditStatusScreen();
                else JOptionPane.showMessageDialog(mFrame, "Не выбрано время в расписании.");
            }
            if(e.getSource() == mStatusCancelButton)
                mStatusEditFrame.setVisible(false);
            if(e.getSource() == mStatusOkButton) {
                if(mStatusEditList.getSelectedValue() == null)
                    JOptionPane.showMessageDialog(mStatusEditFrame, "Не выбран статус заказа");
                else {
                    mStatusEditFrame.setVisible(false);
                    mInfoLabel.setText("");
                    String result = mController.editStatus(mSheduleTable.getModel().getValueAt(mSheduleTable.getSelectedRow(), mSheduleTable.getSelectedColumn()) +
                            "-" + mSheduleTable.getColumnName(mSheduleTable.getSelectedColumn()), mStatusEditList.getSelectedValue().toString(), mStatusDescriptionText.getText());
                    JOptionPane.showMessageDialog(mFrame, result);
                }
            }
            if(e.getSource() == mDeleteButton) {
                if(mSheduleTable.getSelectedColumn() == -1 || mSheduleTable.getSelectedRow() == -1)
                    JOptionPane.showMessageDialog(mStatusEditFrame, "Не выбрано время в расписании");
                else {
                    int confirm = JOptionPane.showConfirmDialog(mFrame, "Вы уверены?", "", JOptionPane.YES_NO_OPTION);
                    if (confirm == JOptionPane.YES_OPTION) {
                        String result = mController.deleteRecord(mSheduleTable.getModel().getValueAt(mSheduleTable.getSelectedRow(), mSheduleTable.getSelectedColumn())
                                + "-" + mSheduleTable.getColumnName(mSheduleTable.getSelectedColumn()));
                        fillTable();
                        JOptionPane.showMessageDialog(mFrame, result);
                    }
                }
            }
            if (e.getSource() == mChangeTimeButton) {
                if (!isTimeChanging) {
                    if (mSheduleTable.getSelectedColumn() == -1 || mSheduleTable.getSelectedRow() == -1) {
                        JOptionPane.showMessageDialog(mStatusEditFrame, "Не выбрано время в расписании.");
                    } else {
                        isTimeChanging = true;
                        mChangeTimeButton.setText("Отменить смену времени");
                        oldTime = mSheduleTable.getModel().getValueAt(mSheduleTable.getSelectedRow(), mSheduleTable.getSelectedColumn())
                                + "-" + mSheduleTable.getColumnName(mSheduleTable.getSelectedColumn());
                    }
                } else {
                    isTimeChanging = false;
                    mChangeTimeButton.setText("Сменить время");
                }
            }
            if (e.getSource() == mRefreshButton) {
                fillTable();
            }
        }
    }
}
