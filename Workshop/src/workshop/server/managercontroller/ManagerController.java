package workshop.server.managercontroller;

import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;
import workshop.server.common.DataAccess;
import workshop.server.common.IDataAccess;

public class ManagerController implements IManagerController {
    
    IDataAccess mDataAccess = DataAccess.getInstance();
    
    public ManagerController() {
    }
    
    @Override
    public void processRequest(Request request, ISocketConnection socketConnection) {
        switch (request.requestType) {
            case Constants.GET_SHEDULE: {
                getShedule(socketConnection);
                break;
            }
            case Constants.GET_INFO: {
                getInfo(socketConnection, request.body);
                break;
            }
            case Constants.SET_STATUS: {
                setStatus(socketConnection, request.body.substring(0, request.body.indexOf(';')),
                        request.body.substring(request.body.indexOf(';') + 1, request.body.indexOf('+')),
                        request.body.substring(request.body.indexOf('+') + 1));
                break;
            }
            case Constants.DELETE_RECORD: {
                deleteRecord(socketConnection, request.body);
                break;
            }
            case Constants.CHANGE_TIME: {
                changeTime(socketConnection, request.body.substring(0, request.body.indexOf(';')), request.body.substring(request.body.indexOf(';') + 1));
                break;
            }
            default: {
                break;
            }
        }
    }   
    
    public void sendFail(ISocketConnection socketConnection) {
        Request response = new Request();
        response.senderType = Constants.SERVER;
        response.requestType = Constants.RESPONSE;
        response.body = Constants.FAIL;
        socketConnection.send(response);
    }

    @Override
    public void getShedule(ISocketConnection socketConnection) {
        String times = mDataAccess.getShedule();
        if (times != null) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.GET_SHEDULE;
            response.body = times;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    @Override
    public void getInfo(ISocketConnection socketConnection, String time) {
        String info = mDataAccess.getInfo(time);
        if (info != null) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.GET_INFO;
            response.body = info;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    public void setStatus(ISocketConnection socketConnection, String time, String status, String statusDescription) {
        String result = mDataAccess.setStatus(time, status, statusDescription);
        if (result != null) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.SET_STATUS;
            response.body = result;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    public void deleteRecord(ISocketConnection socketConnection, String time) {
        String result = mDataAccess.deleteRecord(time);
        if (result != null) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.DELETE_RECORD;
            response.body = result;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    public void changeTime(ISocketConnection socketConnection, String oldTime, String newTime) {
        String result = mDataAccess.changeTime(oldTime, newTime);
        if (result != null) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.DELETE_RECORD;
            response.body = result;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }
}
