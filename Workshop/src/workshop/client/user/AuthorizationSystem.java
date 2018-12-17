package workshop.client.user;

import workshop.common.Constants;
import workshop.common.ISocketConnection;
import workshop.common.Request;

public class AuthorizationSystem implements IAuthorizationSystem{
    ISocketConnection mSocketConnection;

    public AuthorizationSystem(final ISocketConnection socketConnection) {
        mSocketConnection = socketConnection;
    }

    @Override
    public boolean authorization(final String username, final String password) {
        Request request = new Request();
        request.senderType = Constants.USER;
        request.requestType = Constants.AUTHORIZATION;
        request.body = username + ":" + password;
        Request response = mSocketConnection.sendAndGetResponse(request);

        return response.senderType.equals(Constants.SERVER) &&
                response.requestType.equals(Constants.AUTHORIZATION) &&
                response.body.equals(Constants.SUCCESS);
    }
    
}
