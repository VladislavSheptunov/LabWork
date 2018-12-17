package workshop.client.head;

public interface IHController {
    /**
     * Gets server log.
     * @return log.
     */
    public String getLog();
    
    public boolean clearLog();
}
