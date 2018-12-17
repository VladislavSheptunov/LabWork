package workshop.common;

import java.io.Closeable;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.InetAddress;
import java.net.Socket;

public class SocketConnection implements ISocketConnection{
    private boolean mConnected = false;
    private Socket mSocket = null;
    private ObjectInputStream mInputStream = null;
    private ObjectOutputStream mOutputStream = null;
    
    public SocketConnection() {
    }

    @Override
    public boolean connect(final String ip, final int port) { 
		try {
            mConnected = false;
            InetAddress ipAddress = InetAddress.getByName(ip);            
            mSocket = new Socket(ipAddress, port);    
            mOutputStream = new ObjectOutputStream(mSocket.getOutputStream());
            mInputStream = new ObjectInputStream(mSocket.getInputStream());
            mConnected = true;
            return true;
        } catch (IOException e) {
            System.out.println(e.toString());
            return false;
        }
    }
    
    @Override
    public boolean connect(final Socket socket) { 
		try {
            mConnected = false;           
            mSocket = socket;         
            mInputStream = new ObjectInputStream(mSocket.getInputStream());
            mOutputStream = new ObjectOutputStream(mSocket.getOutputStream());
            mConnected = true;
            return true;
        } catch (IOException e) {
            System.out.println(e.toString());
            return false;
        }
    }

    @Override
    public void disconnect() {
        close(mOutputStream);
        close(mInputStream);
        close(mSocket); 
        mConnected = false;
    }
    
    private void close(final Closeable closeable) { 
	if (closeable != null) {
            try {
                closeable.close();
            } catch (IOException e) {
                System.out.println(e.toString());
            } 
        }
    }   

    @Override
    public Request sendAndGetResponse(final Request request) { 
		Request response = new Request();
        response.senderType = Constants.SERVER;
        response.requestType = Constants.RESPONSE;
        response.body = Constants.FAIL;
        if (mConnected) {
            try {
                mOutputStream.writeObject(request);
                Request serverResponse = (Request)mInputStream.readObject();
                response = serverResponse;
            } catch (IOException | ClassNotFoundException e) {
                System.out.println(e.toString());
            }
        }   
        return response;
    }

    @Override
    public boolean send(final Request request) { 
		if (mConnected) {
            try {
                mOutputStream.writeObject(request);
                return true;
            } catch (IOException e) {
                System.out.println(e.toString());
            }
        }
       return false;
    }

    @Override
    public Request getResponse() {
        Request response = new Request();
        response.senderType = Constants.SERVER;
        response.requestType = Constants.RESPONSE;
        response.body = Constants.FAIL;
        if (mConnected) {
            try {
                Request serverResponse = (Request)mInputStream.readObject();
                response = serverResponse;
            } catch (IOException | ClassNotFoundException e) {
                response.body = Constants.CONNECTION_REFUSED;
                System.out.println(e.toString());
            }
        }    
        return response;
    }

    @Override
    public boolean isConnected() {
        return mConnected;
    }
}
