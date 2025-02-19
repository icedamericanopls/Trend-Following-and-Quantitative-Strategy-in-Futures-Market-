//import java.util.Date;

public class data_i {
    public String date;
    public String time;
    public Double open;
    public Double high;
    public Double low;
    public Double close;
    public Double volume;


    //constructor
    public data_i(String date, String time, Double open, Double high, Double low, Double close, Double volume){
        this.date = date;
        this.time = time;
        this.open = open;
        this.high = high;
        this.low = low;
        this.close = close;
        this.volume = volume;
    }

    public String toString() {
        return "Data_i{" +
                "date=" + date + '\'' +
                ", time=" + time + '\'' +
                ", open price =" + open + '\'' +
                ", highest price =" + high + '\'' +
                ", lowest price =" + low + '\'' +
                ", close price =" + close + '\'' +
                ", volume =" + volume + '\'' +
                '}';
    }
    
}

    

