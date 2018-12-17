package workshop.server.usercontroller;

import workshop.common.ISocketConnection;
import workshop.common.Request;

public interface IUserController {
    
    /**
     * Process the requeest.
     * @param request is input request.
     * @param socketConection is socket connection with the client.
     */
    public void processRequest(final Request request, final ISocketConnection socketConection);
    
    /**
     * Registration of user with current login and password.
     * @param username is user's username.
     * @param password is user's password.
     * @param socketConnection is socket connection with user.
     */
    public void registration(final String username, final String password,
            ISocketConnection socketConnection);
    
    /**
     * Authorization of user with current login and password.
     * @param username is user's username.
     * @param password is user's password.
     * @param socketConnection is socket connection with user.
     */
    public void authorization(final String username, final String password,
            ISocketConnection socketConnection);
    
    /**
     * Gets shedule.
     * @param socketConnection is socket connection with user.
     */
    public void getShedule(final ISocketConnection socketConnection);
    
    /**
     * Makes order.
     * @param username is name of user.
     * @param description is description of the proplem.
     * @param phone is phine number.
     * @param time is time from shedule.
     * @param socketConnection is socket connection with the client.
     */
    public void makeOrder(final String username, final String description, 
            final String phone, final String time,
            ISocketConnection socketConnection);
    
    /**
     * Gets orders of current user.
     * @param username is user's name.
     * @param socketConnection is socket connection with the client.
     */
    public void getOrders(final String username, final ISocketConnection socketConnection);
    
    /**
     * Sets status to order as accepted.
     * @param username is name of user.
     * @param time is time from shedule.
     * @param socketConnection is socket connection with the client.
     */
    public void acceptOrder(final String username, final String time, 
            final ISocketConnection socketConnection);
    
    /**
     * Sets status to order as rejected.
     * @param username is name of user.
     * @param time is time from shedule.
     * @param socketConnection is socket connection with the client.
     */
    public void rejectOrder(final String username, final String time,
             final ISocketConnection socketConnection);
}
