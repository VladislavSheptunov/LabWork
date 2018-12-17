package workshop.server.common;

import java.net.Socket;
import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;
import workshop.common.SocketConnection;
import workshop.server.headcontroller.HeadController;
import workshop.server.managercontroller.IManagerController;
import workshop.server.managercontroller.ManagerController;
import workshop.server.usercontroller.IUserController;
import workshop.server.usercontroller.UserController;

public class ClientThread extends Thread {
    private final ISocketConnection mSocketConnection;
    private final Socket mClientSocket;
    private boolean mBreaked = false;
    
    private final IUserController mUserController;
    private final IManagerController mManagerController;
    private final HeadController mHeadController;
    
    public ClientThread(final Socket socket) { 
		mClientSocket = socket;
        mSocketConnection = new SocketConnection();
        mUserController = new UserController();
        mManagerController = new ManagerController();
        mHeadController = new HeadController();
    }
    
    @Override
    public void interrupt() {
        mBreaked = true;
    }
    
    @Override
    public void run() {
        if (mSocketConnection.connect(mClientSocket)) {
            while (!mBreaked) {
                Request request = (Request)mSocketConnection.getResponse();
                if (request != null) {
                    if (request.body.equals(Constants.CONNECTION_REFUSED)) {
                        mSocketConnection.disconnect();
                        mBreaked = true;
                    } else {
                        switch (request.senderType) {
                            case Constants.USER: {
                                mUserController.processRequest(request, mSocketConnection);
                                break;
                            }
                            case Constants.MANAGER: {
                                mManagerController.processRequest(request, mSocketConnection);
                                break;
                            }
                            case Constants.HEAD: {
                                mHeadController.processRequest(request, mSocketConnection);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}
