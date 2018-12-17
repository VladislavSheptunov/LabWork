package workshop.server.usercontroller;

import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;
import workshop.server.common.DataAccess;
import workshop.server.common.IDataAccess;

public class UserController implements IUserController {
    
    IDataAccess mDataAccess = DataAccess.getInstance();
    
    public UserController() {
    }
    
    @Override
    public void registration(final String username, final String password,
            final ISocketConnection socketConnection) {
        if (mDataAccess.addUser(username, password)) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.REGISTRATION;
            response.body = Constants.SUCCESS;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    @Override
    public void authorization(final String username, final String password,
            final ISocketConnection socketConnection) {
        if (mDataAccess.checkUser(username, password)) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.AUTHORIZATION;
            response.body = Constants.SUCCESS;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    @Override
    public void processRequest(final Request request, final ISocketConnection socketConnection) {
        switch (request.requestType) {
            case Constants.REGISTRATION: {
                String[] body = request.body.split(":", 2);
                if (body.length == 2) {
                    registration(body[0], body[1], socketConnection);
                } else {
                    sendFail(socketConnection);
                }
                break;
            }
            case Constants.AUTHORIZATION: {
                String[] body = request.body.split(":", 2);
                if (body.length == 2) {
                    authorization(body[0], body[1], socketConnection);
                } else {
                    sendFail(socketConnection);
                }
                break;
            }
            case Constants.GET_SHEDULE: {
                getShedule(socketConnection);
                break;
            }
            case Constants.MAKE_ORDER: {
                String[] body = request.body.split("=", 4);
                if (body.length == 4) {
                    makeOrder(body[0], body[1], body[2], body[3], socketConnection);
                } else {
                    sendFail(socketConnection);
                }
                break;
            }
            case Constants.GET_ORDERS: {
                getOrders(request.body, socketConnection);
                break;
            }
            case Constants.ACCEPT_ORDER: {
                String[] body = request.body.split("=", 2);
                if (body.length == 2) {
                    acceptOrder(body[0], body[1], socketConnection);
                } else {
                    sendFail(socketConnection);
                }
                break;
            }
            case Constants.REJECT_ORDER: {
                String[] body = request.body.split("=", 2);
                if (body.length == 2) {
                    rejectOrder(body[0], body[1], socketConnection);
                } else {
                    sendFail(socketConnection);
                }
                break;
            }
        }
    }   
    
    private void sendFail(final ISocketConnection socketConnection) {
        Request response = new Request();
        response.senderType = Constants.SERVER;
        response.requestType = Constants.RESPONSE;
        response.body = Constants.FAIL;
        socketConnection.send(response);
    }

    @Override
    public void getShedule(final ISocketConnection socketConnection) {
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
    public void makeOrder(final String username, final String description, 
            final String phone, final String time,
            ISocketConnection socketConnection) {
        if (mDataAccess.makeOrder(username, description, phone, time)) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.MAKE_ORDER;
            response.body = Constants.SUCCESS;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    @Override
    public void getOrders(final String username, final ISocketConnection socketConnection) {
        String orders = mDataAccess.getOrders(username);
        if (orders != null) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.GET_ORDERS;
            response.body = orders;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    @Override
    public void acceptOrder(final String username, final String time, 
            final ISocketConnection socketConnection) {
        if (mDataAccess.acceptOrder(username, time)) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.ACCEPT_ORDER;
            response.body = Constants.SUCCESS;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }

    @Override
    public void rejectOrder(final String username, final String time, 
             final ISocketConnection socketConnection) {
        if (mDataAccess.rejectOrder(username, time)) {
            Request response = new Request();
            response.senderType = Constants.SERVER;
            response.requestType = Constants.REJECT_ORDER;
            response.body = Constants.SUCCESS;
            socketConnection.send(response);
        } else {
            sendFail(socketConnection);
        }
    }
}
