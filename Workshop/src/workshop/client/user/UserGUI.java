package workshop.client.user;

import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.LinkedList;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.WindowConstants;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import workshop.client.common.IGUI;
import workshop.common.Constants;
import workshop.client.common.SheduleTableRenderer;

public class UserGUI implements IGUI {
    private static final int WIDTH = 500;
    private static final int HEIGHT = 500;
    
    private static final String TITLE = "Клиент пользователя";
    private static final String REGISTRATION_TITLE = "Регистрация";
    private static final String LOGIN_TITLE = "Войти";
    private static final String LOGOUT_TITLE = "Выйти";
    private static final String ACCEPT_TITLE = "Принять";
    private static final String REJECT_TITLE = "Отклонить";
    private static final String CANCEL_TITLE = "Назад";
    private static final String REFRESH_TITLE = "Обновить расписание";
    private static final String BACK_TITLE = "Назад";
    private static final String PHONE_TITLE = "Телефон";
    private static final String DESCRIPTION_TITLE = "Описание заказа";
    private static final String STATUS_DESCRIPTION_TITLE = "Описание статуса";
    private static final String SHEDULE_TITLE = "Расписание";
    private static final String USERNAME_TITLE = "Имя";
    private static final String PASSWORD_TITLE = "Пароль";
    private static final String SHOW_TITLE = "Показать мои заказы";
    private static final String SIGN_UP_TITLE = "Записаться";
    
    private final JFrame mFrame;    
    private final JLabel mUsernameLabel;
    private final JLabel mPasswordLabel;
    private final JTextField mUsernameField;
    private final JTextField mPasswordField;
    private final JButton mRegistrationButton;
    private final JButton mLoginButton;
    private final JButton mShowStatusButton;
    private final JButton mSignUpButton;   
    private final JButton mLogoutButton; 
    private final JButton mAcceptButton;
    private final JButton mCancelButton;
    private final JButton mBackButton;
    private final JButton mRefreshButton;
    private final JButton mOrderAcceptButton;
    private final JButton mOrderRejectButton;
    private final JTable mSheduleTable;
    private final JTable mOrderTable;
    private final JTextArea mDescription;
    private final JTextArea mOrderDescription;
    private final JTextArea mStatusDescription;
    private final JLabel mDescriptionLabel;
    private final JLabel mStatusDescriptionLabel;
    private final JLabel mSheduleLabel;
    private final JLabel mPhoneLabel;
    private final JLabel mOrderDescriptionLabel;
    private final JTextField mPhone;
    private final JScrollPane mScrollPane;
    private final JScrollPane mOrderScrollPane;
    
    private final ClickListener mClickListener;
    
    private final IUController mController;
    
    private String mUserName = null;
    private Point[] mBusyTimes = null;
    private final LinkedList<Pair> mDescriptions = new LinkedList();
    
    public UserGUI(final IUController controller) {
        mController = controller;

        mFrame = new JFrame(TITLE);
        mFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        mFrame.setSize(WIDTH, HEIGHT);
        mFrame.setLayout(null);
        
        mUsernameLabel = new JLabel(USERNAME_TITLE);
        mUsernameLabel.setBounds(50, 50, 80, 20);
        
        mPasswordLabel = new JLabel(PASSWORD_TITLE);
        mPasswordLabel.setBounds(50, 90, 80, 20);
        
        mUsernameField = new JTextField();
        mUsernameField.setBounds(150, 50, 100, 20);
        
        mPasswordField = new JTextField();
        mPasswordField.setBounds(150, 90, 100, 20);
        
        mClickListener = new ClickListener();
        
        mLoginButton = new JButton(LOGIN_TITLE);
        mLoginButton.setBounds(50, 130, 200, 30);
        mLoginButton.addActionListener(mClickListener);
        
        mRegistrationButton = new JButton(REGISTRATION_TITLE);
        mRegistrationButton.setBounds(50, 170, 200, 30);
        mRegistrationButton.addActionListener(mClickListener);
        
        mShowStatusButton = new JButton(SHOW_TITLE);
        mShowStatusButton.setBounds(50, 50, 200, 30);
        mShowStatusButton.addActionListener(mClickListener);
        
        mSignUpButton = new JButton(SIGN_UP_TITLE);
        mSignUpButton.setBounds(50, 90, 200, 30);
        mSignUpButton.addActionListener(mClickListener);
        
        mLogoutButton = new JButton(LOGOUT_TITLE);
        mLogoutButton.setBounds(50, 130, 200, 30);
        mLogoutButton.addActionListener(mClickListener);
        
        mAcceptButton = new JButton(ACCEPT_TITLE);
        mAcceptButton.setBounds(50, 20, 200, 30);
        mAcceptButton.addActionListener(mClickListener);     
        
        mCancelButton = new JButton(CANCEL_TITLE);
        mCancelButton.setBounds(260, 20, 200, 30);
        mCancelButton.addActionListener(mClickListener);
        
        mOrderAcceptButton = new JButton(ACCEPT_TITLE);
        mOrderAcceptButton.setBounds(50, 380, 200, 30);
        mOrderAcceptButton.addActionListener(mClickListener);
        
        mOrderRejectButton = new JButton(REJECT_TITLE);
        mOrderRejectButton.setBounds(260, 380, 200, 30);
        mOrderRejectButton.addActionListener(mClickListener);
        
        mBackButton = new JButton(BACK_TITLE);
        mBackButton.setBounds(50, 20, 200, 30);
        mBackButton.addActionListener(mClickListener);
        
        mRefreshButton = new JButton(REFRESH_TITLE);
        mRefreshButton.setBounds(250, 20, 200, 30);
        mRefreshButton.addActionListener(mClickListener);
        
        mDescriptionLabel = new JLabel(DESCRIPTION_TITLE);
        mDescriptionLabel.setBounds(50, 60, 150, 20);
        
        mStatusDescriptionLabel = new JLabel(STATUS_DESCRIPTION_TITLE);
        mStatusDescriptionLabel.setBounds(50, 290, 200, 30);
        
        mOrderDescriptionLabel = new JLabel(DESCRIPTION_TITLE);
        mOrderDescriptionLabel.setBounds(50, 210, 150, 20);
        
        mDescription = new JTextArea();
        mDescription.setBounds(50, 80, 400, 50);
        mDescription.setLineWrap(true);
        
        mStatusDescription = new JTextArea();
        mStatusDescription.setBounds(50, 320, 400, 50);
        
        mOrderDescription = new JTextArea();
        mOrderDescription.setBounds(50, 235, 400, 50);
        mOrderDescription.setLineWrap(true);
        
        mSheduleLabel = new JLabel(SHEDULE_TITLE);
        mSheduleLabel.setBounds(50, 190, 80, 20);
        
        mSheduleTable = new JTable();
        mSheduleTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        mSheduleTable.setBounds(50, 210, 400, 150);
        
        mOrderTable = new JTable();
        mOrderTable.setBounds(50, 50, 400, 150);
        mOrderTable.addMouseListener(new MouseAdapter(){
            @Override
            public void mouseClicked(MouseEvent e) {
                int row = mOrderTable.rowAtPoint(e.getPoint());
                if (!mDescriptions.get(row).first.equals("null")) {
                    mOrderDescription.setText(mDescriptions.get(row).first);
                } else {
                    mOrderDescription.setText("");
                }
                if (!mDescriptions.get(row).second.equals("null")) {
                    mStatusDescription.setText(mDescriptions.get(row).second);
                } else {
                    mStatusDescription.setText("");
                }
                if (mOrderTable.getValueAt(row, 2).equals(Constants.STATUS_NEED_FOR_DETAILS) ||
                        mOrderTable.getValueAt(row, 2).equals(Constants.STATUS_USER_AGREE) ||
                        mOrderTable.getValueAt(row, 2).equals(Constants.STATUS_USER_DENIED)) {
                    mOrderAcceptButton.setEnabled(true);
                    mOrderRejectButton.setEnabled(true);
                } else {
                    mOrderAcceptButton.setEnabled(false);
                    mOrderRejectButton.setEnabled(false);
                }
            }
        });
        
        mPhone = new JTextField();
        mPhone.setBounds(150, 140, 100, 20);
        
        mPhoneLabel = new JLabel(PHONE_TITLE);
        mPhoneLabel.setBounds(50, 140, 100, 20);
        
        mScrollPane = new JScrollPane(mSheduleTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, 
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        mScrollPane.setBounds(50, 210, 400, 150);
        
        mOrderScrollPane = new JScrollPane(mOrderTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, 
                JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        mOrderScrollPane.setBounds(50, 50, 400, 150);      
    }
    
    @Override
    public void show() {
        mFrame.setVisible(true);
        showStartScreen();
    } 
    
    public void showStartScreen() {
        clear();
        mFrame.add(mUsernameLabel);
        mFrame.add(mPasswordLabel);
        mFrame.add(mUsernameField);
        mFrame.add(mPasswordField);
        mFrame.add(mRegistrationButton);
        mFrame.add(mLoginButton);
        mFrame.repaint();
    }
    
    public void showAuthorizationScreen() {
        clear();
        mFrame.add(mSignUpButton);
        mFrame.add(mShowStatusButton);
        mFrame.add(mLogoutButton);
        mFrame.repaint();
    }
    
    public void showSignUpScreen() {
        clear();
        mFrame.add(mAcceptButton);
        mFrame.add(mCancelButton);
        mFrame.add(mScrollPane);
        mFrame.add(mDescription);
        mFrame.add(mDescriptionLabel);
        mFrame.add(mSheduleLabel);
        mFrame.add(mPhone);
        mFrame.add(mPhoneLabel);
        fillTable();
        mFrame.repaint();
    }
    
    public void showOrderScreen() {
        clear();
        mOrderDescription.setText("");
        mStatusDescription.setText("");
        mFrame.add(mOrderScrollPane);
        mFrame.add(mBackButton);
        mFrame.add(mOrderDescriptionLabel);
        mFrame.add(mOrderDescription);
        mFrame.add(mStatusDescriptionLabel);
        mFrame.add(mStatusDescription);
        mFrame.add(mOrderAcceptButton);
        mFrame.add(mOrderRejectButton);
        mFrame.add(mRefreshButton);
        mOrderAcceptButton.setEnabled(false);
        mOrderRejectButton.setEnabled(false);
        fillOrderTable();
        mFrame.repaint();
    }
    
    public void clear() {
        mFrame.remove(mUsernameLabel);
        mFrame.remove(mPasswordLabel);
        mFrame.remove(mUsernameField);
        mFrame.remove(mPasswordField);
        mFrame.remove(mRegistrationButton);
        mFrame.remove(mLoginButton);
        mFrame.remove(mSignUpButton);
        mFrame.remove(mShowStatusButton);
        mFrame.remove(mLogoutButton);
        mFrame.remove(mCancelButton);
        mFrame.remove(mAcceptButton);
        mFrame.remove(mSheduleTable);
        mFrame.remove(mDescription);
        mFrame.remove(mDescriptionLabel);
        mFrame.remove(mSheduleLabel);
        mFrame.remove(mPhone);
        mFrame.remove(mPhoneLabel);
        mFrame.remove(mScrollPane);
        mFrame.remove(mOrderScrollPane);
        mFrame.remove(mOrderTable);
        mFrame.remove(mBackButton);
        mFrame.remove(mOrderDescriptionLabel);
        mFrame.remove(mOrderDescription);
        mFrame.remove(mStatusDescriptionLabel);
        mFrame.remove(mStatusDescription);
        mFrame.remove(mOrderAcceptButton);
        mFrame.remove(mOrderRejectButton);
        mFrame.remove(mRefreshButton);
        mFrame.repaint();
    }
    
    public void clearTable(JTable table) {
        DefaultTableModel tableModel = (DefaultTableModel) table.getModel();
        tableModel.setRowCount(0);
        tableModel.setColumnCount(0);
    }
    
    public void fillTable() {
        clearTable(mSheduleTable);
        DateFormat dateFormat = new SimpleDateFormat("dd.MM.yyyy");
        Calendar calendar = Calendar.getInstance();
        String currentDate = dateFormat.format(calendar.getTime());
        Object[] columns = new Object[30];
        for (int i = 0; i < 30; i++) {
            columns[i] = (String)currentDate;
            calendar.add(Calendar.DATE, 1);
            currentDate = dateFormat.format(calendar.getTime());
        }
        TableModel tableModel = new DefaultTableModel(columns, 0);
        mSheduleTable.setModel(tableModel);
        DefaultTableModel model = (DefaultTableModel) mSheduleTable.getModel();
        int hour = 10;
        for (int i = 0; i < 12; i++) {
            Object[] row = new Object[30];
            for (int j = 0; j < 30; j++) {
                row[j] = hour + ":00";
            }
            hour++;
            model.addRow(row);
        }
        mSheduleTable.setDefaultEditor(Object.class, null);
        
        String[] shedule = mController.getShedule();
        if (shedule == null) {
            JOptionPane.showMessageDialog(mFrame,
                                "Не удалось получить расписание!");
        } else {
            int size = shedule.length;
            mBusyTimes = new Point[size];
            for (int i = 0; i < size; i++) {
                String[] splited = shedule[i].split("-");
                if (splited.length == 2) {
                    String time = splited[0];
                    String day = splited[1];
                    int colIndex = getColumnIndex(mSheduleTable, day);
                    if (colIndex != -1) {
                        int rowIndex = getRowIndex(time);
                        if (rowIndex != -1) {
                            mBusyTimes[i] = new Point(rowIndex, colIndex);
                        }
                    }
                }
            }
            setTableRenderer(size, mBusyTimes);
        }
    }
    
    private void setTableRenderer(int size, Point[] filledCells) {
        mSheduleTable.setDefaultRenderer(mSheduleTable.getColumnClass(0), 
                new SheduleTableRenderer(size, filledCells));
    }
    
    public void fillOrderTable() {
        clearTable(mOrderTable);
        Object[] columns = new Object[3];
        columns[0] = "Time";
        columns[1] = "Phone";
        columns[2] = "Status";
        
        TableModel tableModel = new DefaultTableModel(columns, 0);
        mOrderTable.setModel(tableModel);
        DefaultTableModel model = (DefaultTableModel) mOrderTable.getModel();
        mSheduleTable.setDefaultEditor(Object.class, null);
        
        String[] orders = mController.getOrders(mUserName);
        if (orders == null) {
            JOptionPane.showMessageDialog(mFrame,
                                "Не удалось получить заказы!");
        } else {
            mDescriptions.clear();
            mStatusDescription.setText("");
            mOrderDescription.setText("");
            int size = orders.length;
            for (int i = 0; i < size; i++) {
                String[] splited = orders[i].split("=");
                if (splited.length == 5) {
                    String time = splited[0];
                    String phone = splited[1];
                    String status = splited[2];
                    String description = splited[3];
                    String statusDescription = splited[4];
                    mDescriptions.add(new Pair(description, statusDescription));
                    Object[] row = new Object[3];
                    row[0] = time;
                    row[1] = phone;
                    row[2] = status;
                    model.addRow(row);
                }
            }
        }
    }
    
    private int getColumnIndex (JTable table, String header) {
        for (int i = 0; i < table.getColumnCount(); i++) {
            if (table.getColumnName(i).equals(header)) return i;
        }
        return -1;
    }
    
    private int getRowIndex (String header) {
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
            if (e.getSource() == mRegistrationButton) {
                if (!mUsernameField.getText().isEmpty()) {
                    if (mPasswordField.getText().length() >= 8) {
                        if (mController.registration(mUsernameField.getText(),
                                mPasswordField.getText())) {
                            JOptionPane.showMessageDialog(mFrame,
                                "Регистрация завершена!");
                        } else {
                            JOptionPane.showMessageDialog(mFrame,
                                "Имя уже занято!");
                        }
                    } else {
                        JOptionPane.showMessageDialog(mFrame,
                                "Пароль должен иметь не менее 8ми символов!");
                    }
                } else {
                    JOptionPane.showMessageDialog(mFrame, "Имя пользователя пустое!");
                }
            } else if (e.getSource() == mLoginButton) {
                if (!mUsernameField.getText().isEmpty()) {
                    if (mPasswordField.getText().length() >= 8) {
                        if (mController.authorization(mUsernameField.getText(),
                                mPasswordField.getText())) {
                            mUserName = mUsernameField.getText();
                            JOptionPane.showMessageDialog(mFrame,
                                "Авторизация завершена!");
                            showAuthorizationScreen();
                        } else {
                            JOptionPane.showMessageDialog(mFrame,
                                "Не верно имя пользователя или пароль!");
                        }
                    } else {
                        JOptionPane.showMessageDialog(mFrame,
                                "Пароль должен иметь не менее 8ми символов!");
                    }
                } else {
                    JOptionPane.showMessageDialog(mFrame, "Имя пользователя пустое!");
                }
            }  else if (e.getSource() == mLogoutButton) {
                showStartScreen();
            } else if (e.getSource() == mSignUpButton) {
                showSignUpScreen();
            } else if (e.getSource() == mCancelButton) {
                showAuthorizationScreen();
            } else if (e.getSource() == mAcceptButton) {
                if (!mDescription.getText().isEmpty()) {
                    if (!mPhone.getText().isEmpty()) {
                        int selectedRow = mSheduleTable.getSelectedRow();
                        int selectedCol = mSheduleTable.getSelectedColumn();
                        if (selectedRow != -1 && selectedCol != -1) {
                            String time = (String)mSheduleTable.getModel()
                                    .getValueAt(selectedRow, selectedCol);
                            if (mBusyTimes != null &&
                                    !Arrays.asList(mBusyTimes).contains(
                                            new Point(selectedRow, selectedCol))) {
                                if (mUserName != null && 
                                        mController.makeOrder(mUserName, mDescription.getText(),
                                                mPhone.getText(), time + "-" +
                                                        mSheduleTable.getColumnName(selectedCol))) {
                                    JOptionPane.showMessageDialog(mFrame,
                                        "Заказ оформлен!");
                                    showSignUpScreen();
                                } else {
                                    JOptionPane.showMessageDialog(mFrame,
                                        "Это время уже заняли!");
                                    fillTable();
                                }
                            } else {
                                JOptionPane.showMessageDialog(mFrame,
                                    "Это время занято!");
                            }
                        } else {
                            JOptionPane.showMessageDialog(mFrame,
                                "Не выбрано время в расписании!");
                        }
                    } else {
                        JOptionPane.showMessageDialog(mFrame,
                                "Телефонный номер не указан!");
                    }
                } else {
                    JOptionPane.showMessageDialog(mFrame, "Описание заказа пустое!");
                }
            } else if (e.getSource() == mShowStatusButton) {
                showOrderScreen();
            } else if (e.getSource() == mBackButton) {
                showAuthorizationScreen();
            } else if (e.getSource() == mOrderAcceptButton) {
                if (mOrderTable.getSelectedRow() >= 0 && mController.acceptOrder(mUserName, 
                        mOrderTable.getModel().getValueAt(mOrderTable.getSelectedRow(), 0).toString())) {
                    fillOrderTable();
                } else {
                    JOptionPane.showMessageDialog(mFrame,
                                        "Вы должны выбрать заказ!");
                }
            } else if (e.getSource() == mOrderRejectButton) {
                if (mOrderTable.getSelectedRow() >= 0 && mController.rejectOrder(mUserName, 
                        mOrderTable.getModel().getValueAt(mOrderTable.getSelectedRow(), 0).toString())) {
                    fillOrderTable();
                } else {
                    JOptionPane.showMessageDialog(mFrame,
                                        "Вы должны выбрать заказ!");
                }
            } else if (e.getSource() == mRefreshButton) {
                fillOrderTable();
            }
        }
    }
}

class Pair {
    public String first;
    public String second;
    
    public Pair(String first, String second) {
        this.first = first;
        this.second = second;
    }
}