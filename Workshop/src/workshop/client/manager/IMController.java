package workshop.client.manager;

public interface IMController {
    
    /**
     * Gets shedule.
     * @return list of times when workshop is busy.
     */
    public String[] getShedule();

    public String getInfo(String time);

    public String editStatus(String time, String newStatus, String statusDescription);

    public String deleteRecord(String time);

    public String changeTime(String oldTime, String newTime);
}
