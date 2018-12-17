package workshop.server.common;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.logging.Level;
import workshop.common.Constants;

public class DataAccess implements IDataAccess {
    private static final IDataAccess instance_ = new DataAccess();
    private Connection mConnection;
    private final ILogger mLogger = Logger.getInstance();
    private final Object mLock = new Object();
    private DataAccess() {
    }
    
    public static IDataAccess getInstance() {
        return instance_;
    }
    
    @Override
    public boolean isConnectionClosed() {
        try {
            return mConnection.isClosed();
        } catch (SQLException ex) {
            java.util.logging.Logger.getLogger(DataAccess.class.getName()).log(Level.SEVERE, null, ex);
        }
        return false;
    }
    
    @Override
    public Connection getConnection() {
        return mConnection;
    }
    
    @Override
    public boolean connect () {
        try{
            mConnection = DriverManager.getConnection(
            "jdbc:sqlite:workshop.db","root", "1234");
            return true;
        }catch(SQLException e){  }
        return false;
    }

    @Override
    public void disconnect() {
        try {
            mConnection.close();
        } catch (SQLException e) {
            System.out.println(e.toString());
        }
    }

    @Override
    public boolean addUser(final String username, final String password) { 
	synchronized (mLock) {
            try {
                String query = "select count(*) from users WHERE username='" + username + "';";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                resultSet.next();
                if (resultSet.getInt(1) == 0) {
                    query = "INSERT INTO users (username, password) \n" +
                   " VALUES ('" + username + "', '" + password + "');";
                    statement.executeUpdate(query);
                    mLogger.log("Добавлен новый пользователь: " + username);
                    return true;
                }
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return false;
        }
    }

    @Override
    public boolean checkUser(final String username, final String password) {
        synchronized (mLock) {
            try {
                String query = "select count(*) from users WHERE username='" + username + "' AND "
                        + "password='" + password + "';";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                resultSet.next();
                if (resultSet.getInt(1) > 0) {
                    return true;
                }
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return false;
        }
    }

    @Override
    public String getShedule() {
        synchronized (mLock) {
            try {
                String result = "";
                String query = "select time from orders;";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                while (resultSet.next()) {
                    result+= resultSet.getString(1) + ";";
                }
                return result;
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return null;
        }
    }

    @Override
    public boolean makeOrder(final String username, final String description, 
            final String phone, final String time) {
		synchronized (mLock) {
            try {
                String query = "select count(*) from orders WHERE time='" + time + "';";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                resultSet.next();
                if (resultSet.getInt(1) > 0) {
                    return false;
                } 
                query = "select id from users WHERE username='" + username + "';";
                resultSet = statement.executeQuery(query);
                resultSet.next();                              
                int uid = resultSet.getInt(1);
 
                query = "INSERT INTO orders (uid, time, description, status, phone) \n" +
                   " VALUES ('" + uid + "', '" + time + "', '" + description + "', '" + 
                        Constants.STATUS_WAIT_FOR_CAR + "', '" + phone + "');";
                statement.executeUpdate(query);
                StringBuilder sb = new StringBuilder();
                sb.append("Был сделан новый заказ пользователем ").append(username).append(".")
                        .append(" Время заказа: ").append(time).append(".")
                        .append(" Телефон: ").append(phone).append(".")
                        .append(" Описание: ").append(description);
                mLogger.log(sb.toString());
                return true;
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return false;
        }
    }

    @Override
    public String getOrders(final String username) {
        synchronized (mLock) {
            try {
                String result = "";               
                String query = "select id from users WHERE username='" + username + "';";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                resultSet.next();                              
                int uid = resultSet.getInt(1);
                query = "select time, phone, status, description,"
                        + " statusdescription from orders where uid='" + uid + "';";
                resultSet = statement.executeQuery(query);
                while (resultSet.next()) {
                    result+= resultSet.getString(1) + "=" + resultSet.getString(2) +
                            "=" + resultSet.getString(3) + "=" + resultSet.getString(4) + 
                            "=" + resultSet.getString(5) +";";
                }
                return result;
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return null;
        }
    }
	
	@Override  
    public boolean acceptOrder(final String username, final String time) {  
        synchronized (mLock) {  
            try {  
                String query = "select id from users WHERE username='" + username + "';";  
                Statement statement = mConnection.createStatement();  
                ResultSet resultSet = statement.executeQuery(query);  
                resultSet.next();                                
                int uid = resultSet.getInt(1);  
                query = "UPDATE orders SET status= '" +   
                        Constants.STATUS_USER_AGREE + "' WHERE uid='" + uid + "' "  
                        + "AND time='" + time + "';";  
                statement.executeUpdate(query);
                StringBuilder sb = new StringBuilder();
                sb.append("Пользователь ").append(username).append(" принял соглашение.")
                        .append("Время заказа: ").append(time).append(".");
                mLogger.log(sb.toString());
                return true;  
            } catch (SQLException e) {  
                System.out.println(e.toString());  
            }  
            return false;  
        }  
    }  
  
    @Override  
    public boolean rejectOrder(final String username, final String time) {  
            synchronized (mLock) {  
            try {                 
                String query = "select id from users WHERE username='" + username + "';";  
                Statement statement = mConnection.createStatement();  
                ResultSet resultSet = statement.executeQuery(query);  
                resultSet.next();                                
                int uid = resultSet.getInt(1);  
                query = "UPDATE orders SET status= '" +   
                        Constants.STATUS_USER_DENIED + "' WHERE uid='" + uid + "' "  
                        + "AND time='" + time + "';";  
                statement.executeUpdate(query);
                StringBuilder sb = new StringBuilder();
                sb.append("Пользователь ").append(username).append(" отклонил соглашение.")
                    .append("Время заказа: ").append(time).append(".");
                mLogger.log(sb.toString());
                return true;  
            } catch (SQLException e) {  
                System.out.println(e.toString());  
            }  
            return false;  
        }  
    } 
 
    @Override
    public String getInfo(final String time) {
        synchronized (mLock) {
            try {
                String query = "select * from orders WHERE time='" + time + "';";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                if(resultSet.next()) {
                    String result = "";
                    result += resultSet.getString(4) + ";" + resultSet.getString(5) + ";" + resultSet.getString(6) + ";" + resultSet.getString(7);
                    query = "select username from users WHERE id='" + resultSet.getString(2) + "';";
                    statement = mConnection.createStatement();
                    resultSet = statement.executeQuery(query);
                    resultSet.next();
                    result = resultSet.getString(1) + ";" + result;
                    return result;
                } else {
                    return "Время не занято";
                }
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return "Ошибка";
        }
    }

    @Override
    public String setStatus(String time, String status, String statusDescription) {
        synchronized (mLock) {
            try {
                String query = "select * from orders WHERE time='" + time + "';";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                if(resultSet.next()) {
                    query = "UPDATE orders SET status='" + status + "', statusdescription='" + statusDescription + "' WHERE time='" + time + "';";
                    statement = mConnection.createStatement();
                    statement.executeUpdate(query);
                    StringBuilder sb = new StringBuilder();
                    sb.append("Статус заказа был изменен.")
                        .append(" Время заказа: ").append(time).append(". ")
                            .append(" Описание статуса: ").append(statusDescription);
                    mLogger.log(sb.toString());
                    return "Статус был изменен";
                } else {
                    return "Нет заказов в это время";
                }
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return "Ошибка";
        }
    }

    @Override
    public String deleteRecord(String time) {
        synchronized (mLock) {
            try {
                String query = "select * from orders WHERE time='" + time + "';";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                if(resultSet.next()) {
                    query = "delete from orders WHERE time='" + time + "';";
                    statement = mConnection.createStatement();
                    statement.executeUpdate(query);
                    StringBuilder sb = new StringBuilder();
                    sb.append("Заказ был удален.")
                        .append(" Время заказа: ").append(time).append(". ");
                    mLogger.log(sb.toString());
                    return "Заказ был удален";
                } else {
                    return "Нет заказов в это время";
                }
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return "Ошибка";
        }
    }

    @Override
    public String changeTime(String oldTime, String newTime) {
        synchronized (mLock) {
            try {
                String query = "select * from orders WHERE time='" + oldTime + "';";
                Statement statement = mConnection.createStatement();
                ResultSet resultSet = statement.executeQuery(query);
                if(resultSet.next()) {
                    query = "select * from orders WHERE time='" + newTime + "';";
                    statement = mConnection.createStatement();
                    resultSet = statement.executeQuery(query);
                    if (!resultSet.next()) {
                        query = "UPDATE orders SET time='" + newTime + "' WHERE time='" + oldTime + "';";
                        statement = mConnection.createStatement();
                        statement.executeUpdate(query);
                        StringBuilder sb = new StringBuilder();
                        sb.append("Время заказа было изменено с ")
                            .append(oldTime).append(" на ").append(newTime);
                        mLogger.log(sb.toString());
                        return "Время заказа изменено";
                    } else {
                        return "Это время уже занято";
                    }
                } else
                    return "Нет заказов в это время";
            } catch (SQLException e) {
                System.out.println(e.toString());
            }
            return "Ошибка";
        }
    }
}
