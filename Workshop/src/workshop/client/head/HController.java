package workshop.client.head;

import workshop.client.manager.MController;
import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;

public class HController extends MController implements IHController{

    public HController(ISocketConnection socketConnection) {
        super(socketConnection);
    }

    @Override
    public String getLog() {
        Request request = new Request();
        request.senderType = Constants.HEAD;
        request.requestType = Constants.GET_LOG;
        request.body = "";
        Request response = mSocketConnection.sendAndGetResponse(request);

        return response.body;
    } 

    @Override
    public boolean clearLog() {
        Request request = new Request();
        request.senderType = Constants.HEAD;
        request.requestType = Constants.CLEAR_LOG;
        request.body = "";
        Request response = mSocketConnection.sendAndGetResponse(request);

        return response.senderType.equals(Constants.SERVER) &&
                response.requestType.equals(Constants.CLEAR_LOG) &&
                response.body.equals(Constants.SUCCESS);
    }
}
