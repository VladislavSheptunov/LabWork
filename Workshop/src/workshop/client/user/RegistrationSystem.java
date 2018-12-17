package workshop.client.user;

import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;

public class RegistrationSystem implements IRegistrationSystem{
    
    ISocketConnection mSocketConnection;
    
    public RegistrationSystem(final ISocketConnection socketConnection) {
        mSocketConnection = socketConnection;
    }

    @Override
    public boolean checkUsername(final String username) {
        Request request = new Request();
        request.senderType = Constants.USER;
        request.requestType = Constants.CHECK_USERNAME;
        request.body = username;
        Request response = mSocketConnection.sendAndGetResponse(request);

        return response.senderType.equals(Constants.SERVER) &&
                response.requestType.equals(Constants.CHECK_USERNAME) &&
                response.body.equals(Constants.SUCCESS);
    }

    @Override
    public boolean registration(final String username, final String password) {
        Request request = new Request();
        request.senderType = Constants.USER;
        request.requestType = Constants.REGISTRATION;
        request.body = username + ":" + password;
        Request response = mSocketConnection.sendAndGetResponse(request);

        return response.senderType.equals(Constants.SERVER) &&
                response.requestType.equals(Constants.REGISTRATION) &&
                response.body.equals(Constants.SUCCESS);
    }
    
}
