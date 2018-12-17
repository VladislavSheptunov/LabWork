package workshop.client.user;


public interface IRegistrationSystem {
    /**
     * Checks if login exist in database.
     * @param username is input username.
     * @return true if login exists.
     */
    public boolean checkUsername(final String username);
    
    /**
     * Registration of user with current login and password.
     * @param username is user's login.
     * @param password is user's password.
     * @return true if registration was successfully completed.
     * Returns false if user with current login exists.
     */
    public boolean registration(final String username, final String password);
}
