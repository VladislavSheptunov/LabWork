package workshop.server.headcontroller;

import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;
import workshop.server.common.ILogger;
import workshop.server.common.Logger;
import workshop.server.managercontroller.ManagerController;

public class HeadController extends ManagerController implements IHeadController {
    private final ILogger mLogger = Logger.getInstance();
    
    public HeadController() {      
    }

    @Override
    public void getLog(ISocketConnection socketConnection) {
        Request response = new Request();
        response.senderType = Constants.SERVER;
        response.requestType = Constants.GET_LOG;
        response.body = mLogger.getLog();
        socketConnection.send(response);
    }

    @Override
    public void clearLog(ISocketConnection socketConnection) {
        if (mLogger.clear()) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.CLEAR_LOG;
            response.body = Constants.SUCCESS;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }
    
    @Override
    public void processRequest(Request request, ISocketConnection socketConnection) {
        super.processRequest(request, socketConnection);
        switch (request.requestType) {
            case Constants.GET_LOG: {
                getLog(socketConnection);
                break;
            }
            case Constants.CLEAR_LOG: {
                clearLog(socketConnection);
                break;
            }
        }
    }
}
