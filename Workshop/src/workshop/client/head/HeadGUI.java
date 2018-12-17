package workshop.client.head;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import workshop.client.manager.ManagerGUI;

public class HeadGUI extends ManagerGUI {
    private static final String TITLE = "Клиент начальника";
    private static final String SHOW_LOG_TITLE = "Показать журнал";
    private static final String CLEAR_LOG_TITLE = "Очистить журнал";

    private final JButton mShowLogButton;
    private final JButton mClearLogButton;
    private final JTextArea mLogArea;
    final JScrollPane mScroll;
    private final ClickListener mClickListener = new ClickListener();
    
    private final HController mHController;
    
    public HeadGUI(HController controller) {
        super(controller);
        mHController = controller;
        
        mShowLogButton = new JButton(SHOW_LOG_TITLE);
        mShowLogButton.setBounds(470, 50, 200, 30);
        mShowLogButton.addActionListener(mClickListener);
        
        mClearLogButton = new JButton(CLEAR_LOG_TITLE);
        mClearLogButton.setBounds(670, 50, 200, 30);
        mClearLogButton.addActionListener(mClickListener);
        
        mLogArea = new JTextArea();
        mLogArea.setBounds(470, 85, 400, 150);
        
        mScroll = new JScrollPane(mLogArea, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        mScroll.setBounds(470, 85, 400, 150);
    }
    private class ClickListener implements ActionListener {
        @Override
        public void actionPerformed(ActionEvent e) {
            if (e.getSource() == mShowLogButton) {
                mLogArea.setText(mHController.getLog());
            } else if (e.getSource() == mClearLogButton) {
                mHController.clearLog();
                mLogArea.setText(mHController.getLog());
            } 
        }
    }
    
    @Override
    public void show() {
        super.show();
        mFrame.setTitle(TITLE);
        mFrame.setSize(1000, 500);
        showLogScreeen();
    }
    
    public void showLogScreeen() {
        mFrame.add(mShowLogButton);
        mFrame.add(mClearLogButton);
        mFrame.add(mScroll);
    }
}
