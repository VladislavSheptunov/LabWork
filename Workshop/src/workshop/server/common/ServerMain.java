package workshop.server.common;

public class ServerMain {
    public static void main(String[] args) {
	final IDataAccess dataAccess = DataAccess.getInstance(); 
        if (dataAccess.connect()) {
            final ServerThread serverThread = new ServerThread(); 
            serverThread.start();
        } else {
            System.out.println("Can't connect to the DB");
        }
    } 
}
