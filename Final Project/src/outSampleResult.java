import java.util.ArrayList;

public class outSampleResult {
    ArrayList<String> Date_in;
    ArrayList<String> Date_out;
    ArrayList<String> Time_in;
    ArrayList<String> Time_out;
    ArrayList<Integer> Position_before_exit;
    ArrayList<Double> Price_in;
    ArrayList<Double> Price_out;
    ArrayList<Double> Absolute_PnL;
    ArrayList<Double> Percentage_PnL;
    ArrayList<Double> E;
    ArrayList<Double> DD;
    ArrayList<String> recordOptParameters;
    String startDate;
    String endDate;

    //constructor
    public outSampleResult(ArrayList<String> Date_in, ArrayList<String> Date_out, ArrayList<String> Time_in,  
    ArrayList<String> Time_out, ArrayList<Integer> Position_before_exit, ArrayList<Double> Price_in, 
    ArrayList<Double> Price_out, ArrayList<Double> Absolute_PnL, ArrayList<Double> Percentage_PnL, ArrayList<Double> E,
    ArrayList<Double> DD, ArrayList<String> recordOptParameters, String startDate, String endDate){
        this.Date_in = Date_in;
        this.Date_out = Date_out;
        this.Time_in = Time_in;
        this.Time_out = Time_out;
        this.Position_before_exit = Position_before_exit;
        this.Absolute_PnL = Absolute_PnL;
        this.Percentage_PnL = Percentage_PnL;
        this.Price_in = Price_in;
        this.Price_out = Price_out;
        this.E = E;
        this.recordOptParameters = recordOptParameters;
        this.DD = DD;
        this.startDate = startDate;
        this.endDate = endDate;
    }

}
