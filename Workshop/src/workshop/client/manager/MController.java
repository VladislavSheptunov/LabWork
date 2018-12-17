package workshop.client.manager;

import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;

public class MController implements IMController {

    public final ISocketConnection mSocketConnection;

    public MController(ISocketConnection socketConnection) {
        mSocketConnection = socketConnection;
    }

    @Override
    public String[] getShedule() {
        Request request = new Request();
        request.senderType = Constants.MANAGER;
        request.requestType = Constants.GET_SHEDULE;
        request.body = "";
        Request response = mSocketConnection.sendAndGetResponse(request);

        if (response.body.equals(Constants.FAIL)) {
            return null;
        }
        return response.body.split(";");
    }

    @Override
    public String getInfo(String time) {
        Request request = new Request();
        request.senderType = Constants.MANAGER;
        request.requestType = Constants.GET_INFO;
        request.body = time;
        Request response = mSocketConnection.sendAndGetResponse(request);

        if (response.body.equals(Constants.FAIL)) {
            return null;
        }
        return response.body;
    }

    @Override
    public String editStatus(String time, String newStatus, String statusDescription) {
        Request request = new Request();
        request.senderType = Constants.MANAGER;
        request.requestType = Constants.SET_STATUS;
        request.body = time + ";" + newStatus + "+" + statusDescription;
        Request response = mSocketConnection.sendAndGetResponse(request);
        return response.body;
    }

    @Override
    public String deleteRecord(String time) {
        Request request = new Request();
        request.senderType = Constants.MANAGER;
        request.requestType = Constants.DELETE_RECORD;
        request.body = time;
        Request response = mSocketConnection.sendAndGetResponse(request);
        return response.body;
    }

    @Override
    public String changeTime(String oldTime, String newTime) {
        Request request = new Request();
        request.senderType = Constants.MANAGER;
        request.requestType = Constants.CHANGE_TIME;
        request.body = oldTime + ";" + newTime;
        Request response = mSocketConnection.sendAndGetResponse(request);
        return response.body;
    }
}
