package workshop.server.common;

public interface ILogger {
    /**
     * Adds log message.
     * @param message is message to add to log.
     * @return true if operation was successfully completed.
     */
    public boolean log(String message);

    /**
     * Gets log.
     * @return string of log.
     */
    public String getLog();

    /**
     * Clears log.
     * @return true if operation was successfully completed.
     */
    public boolean clear();
}
