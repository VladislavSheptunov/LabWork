package workshop.client.manager;

import javax.swing.JOptionPane;
import workshop.client.common.IGUI;
import workshop.common.ISocketConnection;
import workshop.common.SocketConnection;

public class MMain {

    private static final String IP = "localhost";
    private static final int PORT = 8888;

    public static void main(String[] args) {
        ISocketConnection socketConnection = new SocketConnection();
        if (socketConnection.connect(IP, PORT)) {
            IMController controller = new MController(socketConnection);
            IGUI managerGUI = new ManagerGUI(controller);
            managerGUI.show();
        } else {
            JOptionPane.showMessageDialog(null,
                    "Нет доступа к серверу!");
        }
    }
}
