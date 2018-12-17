package workshop.common;

import java.io.Serializable;

public class Request implements Serializable {
    public String senderType;
    public String requestType;
    public String body;
}
