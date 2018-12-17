package workshop.server.headcontroller;

import workshop.common.ISocketConnection;

public interface IHeadController {
    /**
     * Gets log.
     * @param socketConection is socket connection with the client.
     */
    public void getLog(ISocketConnection socketConection);

    /**
     * Clears log.
     * @param socketConection is socket connection with the client.
     */
    public void clearLog(ISocketConnection socketConection);
}
