package workshop.client.user;

public interface IAuthorizationSystem {
    /**
     * Authorization of user with current login and password.
     * @param username is user's login.
     * @param password is user's password.
     * @return true if authorization was successfully completed.
     */
    public boolean authorization(final String username, String password);
}
