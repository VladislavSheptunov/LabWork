package workshop.client.user;

import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;

class OrderSystem implements IOrderSystem {
    ISocketConnection mSocketConnection;
    
    public OrderSystem(final ISocketConnection socketConnection) {
        mSocketConnection = socketConnection;
    }

    @Override
    public String[] getShedule() {
        Request request = new Request();
        request.senderType = Constants.USER;
        request.requestType = Constants.GET_SHEDULE;
        request.body = "";
        Request response = mSocketConnection.sendAndGetResponse(request);

        if (response.body.equals(Constants.FAIL)) {
            return null;
        }
        return response.body.split(";");
    }

    @Override
    public boolean makeOrder(final String username, final String description, 
            final String phone, final String time) {
        Request request = new Request();
        request.senderType = Constants.USER;
        request.requestType = Constants.MAKE_ORDER;
        request.body = username + "=" + description + "=" + phone + "=" + time;
        Request response = mSocketConnection.sendAndGetResponse(request);

        return response.senderType.equals(Constants.SERVER) &&
                response.requestType.equals(Constants.MAKE_ORDER) &&
                response.body.equals(Constants.SUCCESS);
    }

    @Override
    public String[] getOrders(final String username) {
        Request request = new Request();
        request.senderType = Constants.USER;
        request.requestType = Constants.GET_ORDERS;
        request.body = username;
        Request response = mSocketConnection.sendAndGetResponse(request);

        if (response.body.equals(Constants.FAIL)) {
            return null;
        }
        return response.body.split(";");
    }

    @Override
    public boolean acceptOrder(final String username, final String time) {
        Request request = new Request();
        request.senderType = Constants.USER;
        request.requestType = Constants.ACCEPT_ORDER;
        request.body = username + "=" + time;
        Request response = mSocketConnection.sendAndGetResponse(request);

        return response.senderType.equals(Constants.SERVER) &&
                response.requestType.equals(Constants.ACCEPT_ORDER) &&
                response.body.equals(Constants.SUCCESS);
    }

    @Override
    public boolean rejectOrder(final String username, final String time) {
        Request request = new Request();
        request.senderType = Constants.USER;
        request.requestType = Constants.REJECT_ORDER;
        request.body = username + "=" + time;
        Request response = mSocketConnection.sendAndGetResponse(request);

        return response.senderType.equals(Constants.SERVER) &&
                response.requestType.equals(Constants.REJECT_ORDER) &&
                response.body.equals(Constants.SUCCESS);
    }
}
