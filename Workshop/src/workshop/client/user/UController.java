package workshop.client.user;

import workshop.common.ISocketConnection;

public class UController implements IUController {
    private final ISocketConnection mSocketConnection;
    private final IRegistrationSystem mRegistrationSystem;
    private final IAuthorizationSystem mAuthorizationSystem;
    private final IOrderSystem mOrderSystem;

    public UController(final ISocketConnection socketConnection) {
        mSocketConnection = socketConnection;
        mRegistrationSystem = new RegistrationSystem(mSocketConnection);
        mAuthorizationSystem = new AuthorizationSystem(mSocketConnection);
        mOrderSystem = new OrderSystem(mSocketConnection);
    }

    @Override
    public boolean registration(final String username, final String password) {
        return mRegistrationSystem.registration(username, password);
    }

    @Override
    public boolean authorization(final String username, final String password) {
        return mAuthorizationSystem.authorization(username, password);
    } 

    @Override
    public String[] getShedule() {
        return mOrderSystem.getShedule();
    }

    @Override
    public boolean makeOrder(final String username, final String description, 
            final String phone, final String time) {
        return mOrderSystem.makeOrder(username, description, phone, time);
    }

    @Override
    public String[] getOrders(final String username) {
        return mOrderSystem.getOrders(username);
    }

    @Override
    public boolean acceptOrder(final String username, final String time) {
        return mOrderSystem.acceptOrder(username, time);
    }

    @Override
    public boolean rejectOrder(final String username, final String time) {
        return mOrderSystem.rejectOrder(username, time);
    }
}
