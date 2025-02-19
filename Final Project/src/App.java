import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;


public class App {
    public static void main(String[] args) throws Exception {
        //read csv file into arraylist of data_i
        //use your own file path
        String filepath = "/Users/xinyizhang/Desktop/math method for financial pricing/final project/SY-5minHLV.csv";
        //String filepath = "/Users/xinyizhang/Desktop/math method for financial pricing/final project/SOYBN-5minHLV.csv";
        List<data_i> dataList = readCSV(filepath);
        //check data_i correct or not
        //System.out.println(dataList.get(0).toString());

        //Generate ChnLen List and StpPct List under given ranges
        List<Double> ChnLen = new ArrayList<>();
        List<Double> StpPct = new ArrayList<>();
        // Generate ChnLen and StpPct arraylist
        for (Double num = 500.0; num <= 10000; num += 10) {
            ChnLen.add(num);
        }
        for (Double num = 0.005; num <= 0.1; num  += 0.001){
            StpPct.add(num);
        }
        //10034, 125864
        final int apprNumberOfIntervals = 45; //SY:45, SOYBN:72
        int totalLength = dataList.size();
        final int inSampleYear = 4;
        final int outSampleMonth = 3;
        final int finalInSampleStart = Math.max(10034, 10001);//最少需要前10000行数据来找optimal ChnLen, SY starts from 10034, SOYBN starts from 10007
        final int finalInSampleEnd = totalLength -1 - outSampleMonth*22*apprNumberOfIntervals - inSampleYear*252*apprNumberOfIntervals; //approximated the last inSampleStart index, need to be changed 
        //Math.min(finalInSampleStart + 10*252*45, totalLength - 1);
        int inSampleStart; 
        int inSampleEnd;
        int outSampleStart; 
        int outSampleEnd; 

        //some lists to store out of sample result
        ArrayList<String> outSampleDateIn = new ArrayList<>();
        ArrayList<String> outSampleDateOut = new ArrayList<>();
        ArrayList<String> outSampleTimeIn = new ArrayList<>();
        ArrayList<String> outSampleTimeOut = new ArrayList<>();
        ArrayList<Double> outSamplePriceIn = new ArrayList<>();
        ArrayList<Double> outSamplePriceOut = new ArrayList<>();
        ArrayList<Double> outSampleAbsolutePnL = new ArrayList<>();
        ArrayList<Double> outSamplePercentPnL = new ArrayList<>();
        ArrayList<Integer> outSamplePositionBeforeExit = new ArrayList<>();
        ArrayList<Double> Equity = new ArrayList<>(); 
        ArrayList<String> outSampleParameters = new ArrayList<>();
        ArrayList<String> statistics = new ArrayList<>();
        statistics.add("date\taverageReturn\tdeviation\tsharpeRatio\tnumberOfTrades\tnumberOfProfitTrades\twinnerPercent\tmaxDrawDown\tnetProfitToMaxDD\tnetEquity\tnetProfit");

        Double E0 = 200000.0;
        int PV = 50;
        int slpg = 65;
        //iterate through each in-sample period
        for (int i = finalInSampleStart; i <= finalInSampleEnd; i += 3*22*apprNumberOfIntervals){ 
            inSampleStart = i;
            inSampleEnd = inSampleStart + inSampleYear*252*apprNumberOfIntervals;
            outSampleStart = inSampleEnd + 1;
            outSampleEnd = Math.min(inSampleEnd + outSampleMonth*22*apprNumberOfIntervals, totalLength - 1);//avoid IndexOutOfBounds

            Double[] optParameters = {0.0, 0.0};
            optParameters = findOptimal2(slpg, PV, E0, dataList, ChnLen, StpPct, inSampleStart, inSampleEnd); //find optimal parameters from in-sample data
            System.out.println(optParameters[0] + " " + optParameters[1]);
            outSampleResult result = getOutSampleResult(slpg, PV, E0, dataList, optParameters, outSampleStart, outSampleEnd);//return out-of-sample data

            //get out of sample satistics
            statistics.add(calculateStatistics(result.startDate, result.endDate, result.Absolute_PnL, result.DD, result.E));
            //add out of sample data into lists
            outSampleDateIn.addAll(result.Date_in);
            outSampleDateOut.addAll(result.Date_out);
            outSampleTimeIn.addAll(result.Time_in);
            outSampleTimeOut.addAll(result.Time_out);
            outSamplePriceIn.addAll(result.Price_in);
            outSamplePriceOut.addAll(result.Price_out);
            outSampleAbsolutePnL.addAll(result.Absolute_PnL);
            outSamplePercentPnL.addAll(result.Percentage_PnL);
            outSamplePositionBeforeExit.addAll(result.Position_before_exit);
            Equity.addAll(result.E);
            outSampleParameters.addAll(result.recordOptParameters);
            //update E0 after each trading period, next trading period will start from E0
            E0 = result.E.get(result.E.size() - 1); 
        }
        
        //print out trade-by-trade table in terminal
        ArrayList<String> writeOutList = new ArrayList<>();
        String title = "Date_In\tTime_In\tDate_Out\tTime_Out\tPosition_Before_Exit\tPrice_In\tPrice_Out\tAbsolute_PnL\tPercent_PnL";
        writeOutList.add(title);
        System.out.println("Date In\t\tTime In\t\tDate Out\tTime Out\tPosition Before Exit\tPrice In\tPrice Out\tAbsolute PnL\tPercent PnL");
        for (int i = 0; i < outSampleDateIn.size(); i++) {
            String temp = String.format("%s\t%s\t%s\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f", outSampleDateIn.get(i), 
            outSampleTimeIn.get(i), outSampleDateOut.get(i), outSampleTimeOut.get(i), outSamplePositionBeforeExit.get(i), outSamplePriceIn.get(i),
            outSamplePriceOut.get(i), outSampleAbsolutePnL.get(i), outSamplePercentPnL.get(i));
            writeOutList.add(temp);
            System.out.printf("%-10s\t%-10s\t%-10s\t%-10s\t%-20d\t$%-10.2f\t$%-10.2f\t$%-10.2f\t%-10.5f\n", outSampleDateIn.get(i), 
            outSampleTimeIn.get(i), outSampleDateOut.get(i), outSampleTimeOut.get(i), outSamplePositionBeforeExit.get(i), outSamplePriceIn.get(i),
            outSamplePriceOut.get(i), outSampleAbsolutePnL.get(i), outSamplePercentPnL.get(i));
        }

        //write into csv file
        //use your own file address
        String outputFilePath = "/Users/xinyizhang/Desktop/math method for financial pricing/final project/outputTable.csv";
        String outputFilePath2 = "/Users/xinyizhang/Desktop/math method for financial pricing/final project/equity.csv";
        String outputFilePath3 = "/Users/xinyizhang/Desktop/math method for financial pricing/final project/outSampleParameters.csv";
        String outputFilePath4 = "/Users/xinyizhang/Desktop/math method for financial pricing/final project/outSampleStatistics.csv";
        //write trade-by-trade into csv
        writeIntoCsv(outputFilePath, writeOutList);
        //convert Double into String
        ArrayList<String> EquityString = new ArrayList<>();
        for (Double d : Equity){
            EquityString.add(String.valueOf(d));
        }
        //write 5-min Equity into csv
        writeIntoCsv(outputFilePath2, EquityString);
        //wrtie optimal parameters into csv
        writeIntoCsv(outputFilePath3, outSampleParameters);
        //write statistics into csv
        writeIntoCsv(outputFilePath4, statistics);
    }

    /**写入csv文档的function
     * filePath: csv文档地址
     * writeOutList: 写入文档的内容
    */
    public static void writeIntoCsv(String filePath, ArrayList<String> writeOutList){
        String outputFilePath = filePath;
        // 指定字符编码
        Charset charset = StandardCharsets.UTF_8;
        // 指定缓存
        int bufferSize = 5 * 1024 * 1024;
        try (BufferedWriter writer = new BufferedWriter(
            new OutputStreamWriter(new FileOutputStream(outputFilePath), charset), bufferSize
        )) {
            for (String datum : writeOutList) {
                writer.write(datum);
                writer.newLine();
                }       
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    
    public static String calculateStatistics(String startDate, String endDate, ArrayList<Double> Absolute_PnL, ArrayList<Double> DD, ArrayList<Double> E){
        //calculate average rate of return
        Double averageReturn = calculateMean(Absolute_PnL);
        Double deviation = Math.sqrt(calculateVariance(Absolute_PnL, averageReturn));
        Double sharpeRatio = averageReturn/deviation;
        Double numberOfTrades = (double) Absolute_PnL.size();
        Double netEquity = E.get(E.size() - 1);
        Double netProfit = sum(Absolute_PnL);
        Double maxDrawDown = getMin_double(DD);
        Double numberOfProfitTrades = (double) countPositiveEntries(Absolute_PnL);
        Double winnerPercent = numberOfProfitTrades/numberOfTrades;
        Double netProfitToMaxDD = netProfit/maxDrawDown;
        String period = startDate + "-" + endDate;
        String statistics = String.format("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", period, averageReturn, deviation,
        sharpeRatio, numberOfTrades, numberOfProfitTrades, winnerPercent, maxDrawDown, netProfitToMaxDD, netEquity, netProfit);

        return statistics;
    }

    public static Double sum(ArrayList<Double> list){
        Double sum = 0.0;
        for (Double num : list){
            sum += num;
        }
        return sum;
    }
    public static int countPositiveEntries(ArrayList<Double> list) {
        int count = 0;
        for (Double num : list) {
            if (num > 0) {
                count++;
            }
        }
        return count;
    }
    public static double calculateMean(ArrayList<Double> list) {
        if (list.isEmpty()) {
            return 0; // Return 0 for an empty list or handle it as needed
        }
        Double sum = 0.0;
        for (Double num : list) {
            sum += num;
        }
        return (double) sum / list.size();
    }

    public static double calculateVariance(ArrayList<Double> list, double mean) {
        if (list.isEmpty()) {
            return 0; // Return 0 for an empty list or handle it as needed
        }
        Double sumSquaredDiffs = 0.0;
        for (Double num : list) {
            Double diff = num - mean;
            sumSquaredDiffs += diff * diff;
        }
        return sumSquaredDiffs / list.size();
    }
    

    /*
     * 找到最佳parameters的function
     */
    public static Double[] findOptimal(int slpg, int PV, Double E0, List<data_i> dataList, List<Double> ChnLen, List<Double> StpPct, int inSampleStart, int inSampleEnd){

        //List<data_i> inSampleData = dataList.subList(inSampleStart, inSampleEnd);
        //List<Double[]> PnLTable = new ArrayList<>();
        Double[] optParameters = {0.0, 0.0};
        Double optRatio = Double.MIN_VALUE;
        int position = 0;
        //Double E0 = 200000.0;
        Hashtable<int[], Double> inSampleResult = new Hashtable<>();

        //iterate through each ChnLen
        for (int i = 0; i < ChnLen.size(); i += 1){
            Double L = ChnLen.get(i);
            System.out.println("calculating length = " + L);

            List<Double> HH = new ArrayList<>();
            List<Double> LL = new ArrayList<>();
            //get HH and LL for given length
            for (int k = inSampleStart; k <= inSampleEnd; k += 1){
                HH.add(getMax(dataList.subList(k - L.intValue(), k-1)));
                LL.add(getMin(dataList.subList(k - L.intValue(), k-1)));
            }

            //iterate through each StpPct
            for (int j = 0; j < StpPct.size(); j += 1){
                Double S = StpPct.get(j);
                System.out.println("calculating PctStp = " + S);
                
                //set initial conditions
                position = 0;
                List<Double> E = new ArrayList<>();
                E.add(E0);
                List<Double> DD = new ArrayList<>();
                DD.add(0.0);
                List<Double> trades = new ArrayList<>();
                trades.add(0.0);
                Double Emax = E0;
                boolean traded, buy, sell, sellShort, buyLong;
                Double delta;
                Double  benchmarkLong = 0.0;
                Double benchmarkShort = 0.0;

                //running through insample and trading
                for (int k= inSampleStart; k <= inSampleEnd; k += 1){
                    traded=false;
                    delta=PV*(dataList.get(k).close-dataList.get(k-1).close)*position;
                    //System.out.println(delta);

                    if (position== 0){
                        //System.out.println(dataList.get(k).high);
                        //System.out.println(HH.get(k - inSampleStart - 1));
                        buy=dataList.get(k).high>=HH.get(k - inSampleStart);
                        sell=dataList.get(k).low<=LL.get(k - inSampleStart);
                        
                        if (buy && sell){
                            delta = -slpg+PV*(LL.get(k - inSampleStart)-HH.get(k- inSampleStart));
                            trades.add(1.0);
                        }
                        else{
                            if(buy){
                                delta = -slpg/2 + PV*(dataList.get(k).close-HH.get(k- inSampleStart ));
                                position= 1;
                                traded=true;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(0.5);
                            }
                
                            if(sell){
                                //System.out.println("k: " + k);
                                //System.out.println("sell :" + sell);
                                delta = -slpg/2 - PV*(dataList.get(k).close-LL.get(k- inSampleStart));
                                //System.out.println("Delta: "+ delta);
                                position=-1;
                                traded=true;
                                benchmarkShort=dataList.get(k).low;
                                //System.out.println("benchmarkShort: " + benchmarkShort);
                                trades.add( 0.5);
                            }
                        }
                    }
                    
            
            
                    if (position== 1 && !traded){
                        sellShort=dataList.get(k).low<=LL.get(k - inSampleStart);
                        sell=dataList.get(k).low<=(benchmarkLong*(1-S));
                            
                        if(sellShort && sell){
                            //copy of sell short
                            if(sellShort){
                                    delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - inSampleStart));
                                    position=-1;
                                    benchmarkShort=dataList.get(k).low;
                                    trades.add(1.0);
                            }
                        }
                        else{
                            if(sell){
                                delta=delta-slpg/2-PV*(dataList.get(k).close-(benchmarkLong*(1-S))); 
                                position=0;
                                trades.add(0.5);
                            }
                                
                            if(sellShort){
                                delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - inSampleStart)); 
                                position=-1;
                                benchmarkShort=dataList.get(k).low;
                                trades.add(1.0);
                            }
                        }
                            
                        benchmarkLong=Math.max(dataList.get(k).high,benchmarkLong);
                    }
            
                    if (position==-1 && !traded){
                        buyLong=dataList.get(k).high>=HH.get(k - inSampleStart);
                        buy=dataList.get(k).high>=(benchmarkShort*(1+S));
                        
                        if(buyLong && buy){
                            //copy of buyLong
                            if(buyLong){
                                delta=delta-slpg+2*PV*(dataList.get(k).close - HH.get(k - inSampleStart));
                                position=1;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(1.0);
                            }
                        }
                        else{
                            if(buy){
                                //System.out.println("k: " + k);
                                //System.out.println("short position - Buy");
                                delta=delta-slpg/2+PV*(dataList.get(k).close-(benchmarkShort*(1+S)));
                                //System.out.println("Delta: "+ delta);
                                position=0;
                                trades.add(0.5);
                            }
                            
                            if(buyLong){
                                delta=delta-slpg+2*PV*(dataList.get(k).close-HH.get(k - inSampleStart));
                                position=1;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(1.0);
                            }
                        }
                        
                        benchmarkShort=Math.min(dataList.get(k).low,benchmarkShort);
                    }
            
                    if (position== 0 && traded){} //do nothing
                    
                    if (position== 1 && traded){} //do nothing
                    
                    if (position==-1 && traded){} //do nothing
                    
                    //update equity
                    if (k-inSampleStart != 0){
                        E.add(E.get(k-inSampleStart-1)+delta);
                    }
                    else{
                        E.set(0, E0+delta);
                    }
                    //calculate drawdown
                    Emax=Math.max(Emax, E.get(k-inSampleStart));
                    DD.add(E.get(k-inSampleStart)-Emax);
                }
                
                //get max drawdown
                Double dd_max = getMin_double(DD);
                //get net profit
                Double netProfit = E.get(E.size() - 1) - E0;
                //System.out.println("InSample Net Profit: " + netProfit);
                //System.out.println("InSample Max DrawDown: " + dd_max);
                //compute ratio
                Double netProfit_ddMax_ratio = netProfit/(-dd_max);
                //System.out.println("netProfit_ddMax_ratio: " + netProfit_ddMax_ratio);
                Double[] L_S = {L,S};
                int[] i_j = {i, j};
                //find the best parameters
                if (netProfit_ddMax_ratio > optRatio){
                    optRatio = netProfit_ddMax_ratio;
                    optParameters = L_S;
                }
                //store each ratio with its corresponding parameters            
                inSampleResult.put(i_j, netProfit_ddMax_ratio);
            }
        }
        return optParameters;
    }

    public static outSampleResult getOutSampleResult(int slpg, int PV, Double E0, List<data_i> dataList, Double[] optParameters, int outSampleStart, int outSampleEnd){
        String startDate = dataList.get(outSampleStart).date + " " + dataList.get(outSampleStart).time;
        String endDate = dataList.get(outSampleEnd).date + " " + dataList.get(outSampleEnd).time;
        System.out.println(startDate);
        System.out.println(endDate);
        //initialize variables
        Double L = optParameters[0];
        Double S = optParameters[1];
        ArrayList<String> recordOptParameters = new ArrayList<String>();
        recordOptParameters.add(dataList.get(outSampleStart).date + " " + dataList.get(outSampleStart).time+"\t"+
        dataList.get(outSampleEnd).date + " " + dataList.get(outSampleEnd).time+ String.format("\t%.2f\t%.4f", L, S));
        System.out.printf("Optimal ChnLen: %.2f, Optimal StpPct: %.4f\n", L, S);
        int position = 0;
        ArrayList<Double> E = new ArrayList<>();
        E.add(E0);
        ArrayList<Double> DD = new ArrayList<>();
        DD.add(0.0);
        ArrayList<Double> trades = new ArrayList<>();
        trades.add(0.0);
        Double Emax = E0;
        boolean traded, buy, sell, sellShort, buyLong;
        Double delta;
        Double  benchmarkLong = 0.0;
        Double benchmarkShort = 0.0;
        Double delta_sum = 0.0;

        //output arraries for trade-by-trade table
        ArrayList<String> Date_in = new ArrayList<>();
        ArrayList<String> Date_out = new ArrayList<>();
        ArrayList<String> Time_in = new ArrayList<>();
        ArrayList<String> Time_out = new ArrayList<>();
        ArrayList<Integer> Position_before_exit = new ArrayList<>();
        ArrayList<Double> Price_in = new ArrayList<>();
        ArrayList<Double> Price_out = new ArrayList<>();
        ArrayList<Double> Absolute_PnL = new ArrayList<>();
        ArrayList<Double> Percentage_PnL = new ArrayList<>();
        
        List<Double> HH = new ArrayList<>();
        List<Double> LL = new ArrayList<>();
        //get HH and LL for given length
        for (int k = outSampleStart; k <= outSampleEnd; k += 1){
            HH.add(getMax(dataList.subList(k - L.intValue(), k-1)));
            LL.add(getMin(dataList.subList(k - L.intValue(), k-1)));
        }

        //running through out-of-sample and trading with optimal parameters
        for (int k= outSampleStart; k <= outSampleEnd; k += 1){
            traded=false;
            delta=PV*(dataList.get(k).close-dataList.get(k-1).close)*position;
            delta_sum += delta; //delta_sum is used to record equity change when we long or short futures
            

            if (position== 0){
                buy=dataList.get(k).high>=HH.get(k - outSampleStart);
                sell=dataList.get(k).low<=LL.get(k - outSampleStart);
         
                if (buy && sell){//在五分钟内同时买入和卖出，假设先买入再卖出
                    //System.out.println(dataList.get(k).date + dataList.get(k).time);
                    //System.out.println("Buy and Sell");
                    delta_sum -= delta;
                    delta = -slpg+PV*(LL.get(k - outSampleStart)-HH.get(k- outSampleStart));//买入价格为当时的HH，卖出价格为当时的LL
                    delta_sum += delta;
                    trades.add(1.0);
                    //buy and sell in the same 5-min interval

                    //进行交易，将交易数据记录进trade-by-trade table
                    //进场
                    Date_in.add(dataList.get(k).date);
                    Time_in.add(dataList.get(k).time);
                    Price_in.add(HH.get(k - outSampleStart));
                    //出场
                    Date_out.add(dataList.get(k).date);
                    Time_out.add(dataList.get(k).time);
                    Position_before_exit.add(position);
                    Price_out.add(LL.get(k-outSampleStart));
                    Absolute_PnL.add(delta_sum);
                    Percentage_PnL.add(Absolute_PnL.get(Absolute_PnL.size() - 1)/E0);
                    delta_sum = 0.0; //结束一段交易，delta_sum归位0
                }
                else{
                    if(buy){ //long only
                        //System.out.println(dataList.get(k).date + dataList.get(k).time);
                        //System.out.println("buy only position0");
                        delta_sum -= delta;
                        delta = -slpg/2 + PV*(dataList.get(k).close-HH.get(k- outSampleStart));//Long position 买入价格为当时的HH
                        delta_sum += delta;
                        position= 1;
                        traded=true;
                        benchmarkLong=dataList.get(k).high;
                        trades.add(0.5);

                        //进行交易，将交易数据记录进trade-by-trade table
                        //进场
                        Date_in.add(dataList.get(k).date);
                        Time_in.add(dataList.get(k).time);
                        Price_in.add(HH.get(k - outSampleStart));
                    }
        
                    if(sell){//short only
                        //System.out.println(dataList.get(k).date + dataList.get(k).time);
                        //System.out.println("Sell ONly position0");
                        delta_sum -= delta;
                        delta = -slpg/2 - PV*(dataList.get(k).close-LL.get(k- outSampleStart));//Short position 卖出价格为当时的LL
                        delta_sum += delta;
                        position= -1;
                        traded=true;
                        benchmarkShort=dataList.get(k).low;
                        trades.add( 0.5);

                        //进行交易，将交易数据记录进trade-by-trade table
                        //进场
                        Date_in.add(dataList.get(k).date);
                        Time_in.add(dataList.get(k).time);
                        Price_in.add(LL.get(k-outSampleStart));
                    }
                }
            }
            
    
    
            if (position== 1 && !traded){//当前为Long position且在此五分钟内没有进行过交易，只能退出平仓或者改为short position，不能继续买入
                sellShort=dataList.get(k).low<=LL.get(k - outSampleStart);
                sell=dataList.get(k).low<=(benchmarkLong*(1-S));
                    
                if(sellShort && sell){
                    if(sellShort){//long positoin变为short position
                        //System.out.println(dataList.get(k).date + dataList.get(k).time);
                        //System.out.println("Sell and SellShort");
                        delta_sum -= delta;
                        delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - outSampleStart));
                        delta_sum += delta;
                        Position_before_exit.add(position);
                        position=-1;
                        benchmarkShort=dataList.get(k).low;
                        trades.add(1.0);

                        //sell holdings
                        //进行交易，将交易数据记录进trade-by-trade table
                        //退出，先卖出持有的futures （结束上一次Long position的交易）
                        Date_out.add(dataList.get(k).date);
                        Time_out.add(dataList.get(k).time);
                        Price_out.add(LL.get(k - outSampleStart)); 
                        //这里感觉退出平仓的价格应该是benchmarkLong*（1- S），但计算delta时使用的是当前的LL， 所以这里也改为LL
                        Absolute_PnL.add(delta_sum);
                        Percentage_PnL.add(Absolute_PnL.get(Absolute_PnL.size() - 1)/E0);
                        delta_sum = 0.0;

                        //进场，再short，开始新一轮交易
                        Date_in.add(dataList.get(k).date);
                        Time_in.add(dataList.get(k).time);
                        Price_in.add(LL.get(k - outSampleStart));//short position 卖出价格为当前的LL
                    }
                }
                else{
                    if(sell){//只平仓，退出交易
                        //System.out.println(dataList.get(k).date+ dataList.get(k).time);
                        //System.out.println("Sell only position1");
                        delta_sum -= delta;
                        delta=delta-slpg/2-PV*(dataList.get(k).close-(benchmarkLong*(1-S)));
                        delta_sum += delta;
                        Position_before_exit.add(position);
                        position=0;
                        trades.add(0.5);

                        //退出
                        Date_out.add(dataList.get(k).date);
                        Time_out.add(dataList.get(k).time);
                        Price_out.add(benchmarkLong*(1-S)); //Long position 卖出价格为benchmarkLong*（1- S）
                        Absolute_PnL.add(delta_sum);
                        Percentage_PnL.add(Absolute_PnL.get(Absolute_PnL.size() - 1)/E0);
                        delta_sum = 0.0;
                    }
                        
                    if(sellShort){//Long position变为Short position，先long position平仓，再short position卖出， 同Sell&SellShort
                        //System.out.println(dataList.get(k).date + dataList.get(k).time);
                        //System.out.println("SellShort");
                        delta_sum -= delta;
                        delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - outSampleStart)); 
                        delta_sum += delta;
                        Position_before_exit.add(position);
                        position=-1;
                        benchmarkShort=dataList.get(k).low;
                        trades.add(1.0);

                        //退出，先卖出持有的futures （结束上一次Long position的交易）
                        Date_out.add(dataList.get(k).date);
                        Time_out.add(dataList.get(k).time);
                        Price_out.add(LL.get(k - outSampleStart));
                        Absolute_PnL.add(delta_sum);
                        Percentage_PnL.add(Absolute_PnL.get(Absolute_PnL.size() - 1)/E0);
                        delta_sum = 0.0;

                        //进场，再short，开始新一轮交易
                        Date_in.add(dataList.get(k).date);
                        Time_in.add(dataList.get(k).time);
                        Price_in.add(LL.get(k - outSampleStart));//short position 卖出价格为当前的LL
                    }
                }
                    
                benchmarkLong=Math.max(dataList.get(k).high,benchmarkLong);//更新benchmarkLong为当前最高价
            }
    
            if (position==-1 && !traded){//当前为Short position且在此五分钟内没有进行过交易，只能退出平仓或者改为Long position，不能继续卖出
                buyLong=dataList.get(k).high>=HH.get(k - outSampleStart);
                buy=dataList.get(k).high>=(benchmarkShort*(1+S));
                
                if(buyLong && buy){
                    if(buyLong){
                        //System.out.println(dataList.get(k).date + dataList.get(k).time);
                        //System.out.println(("Buylong and Buy"));
                        delta_sum -= delta;
                        delta=delta-slpg+2*PV*(dataList.get(k).close - HH.get(k - outSampleStart));
                        delta_sum += delta;
                        Position_before_exit.add(position);
                        position=1;
                        benchmarkLong=dataList.get(k).high;
                        trades.add(1.0);

                        //先退出，买入之前short的futures， 结束short postion的交易
                        Date_out.add(dataList.get(k).date);
                        Time_out.add(dataList.get(k).time);
                        Price_out.add(HH.get(k - outSampleStart));
                        //同样感觉这里的买入价格应该是benchmarkShort*（1+S）， 但计算delta用的是HH，所以这里也输入当时的HH
                        Absolute_PnL.add(delta_sum);
                        Percentage_PnL.add(Absolute_PnL.get(Absolute_PnL.size() - 1)/E0);
                        delta_sum = 0.0; 

                        //再进场，买入，开始long position的交易
                        Date_in.add(dataList.get(k).date);
                        Time_in.add(dataList.get(k).time);
                        Price_in.add(HH.get(k - outSampleStart));//long position买入的价格为当时的HH
                    }

                    
                }
                else{
                    if(buy){//退出交易，只买入short position的futures
                        //System.out.println(dataList.get(k).date + dataList.get(k).time);
                        //System.out.println("Buy positoin-1");
                        delta_sum -= delta;
                        delta=delta-slpg/2+PV*(dataList.get(k).close-(benchmarkShort*(1+S)));
                        delta_sum += delta;
                        Position_before_exit.add(position);
                        position=0;
                        trades.add(0.5);

                        //退出，结束short position
                        Date_out.add(dataList.get(k).date);
                        Time_out.add(dataList.get(k).time);
                        Price_out.add(benchmarkShort*(1+S)); //short position买入价格为benchmarkShort*（1+S）
                        Absolute_PnL.add(delta_sum);
                        Percentage_PnL.add(Absolute_PnL.get(Absolute_PnL.size() - 1)/E0);
                        delta_sum = 0.0;
                    }
                    
                    if(buyLong){ //同buy&buyLong
                        //System.out.println(dataList.get(k).date + dataList.get(k).time);
                        //System.out.println("Buylong");
                        delta_sum -= delta;
                        delta=delta-slpg+2*PV*(dataList.get(k).close-HH.get(k - outSampleStart));
                        delta_sum += delta;
                        Position_before_exit.add(position);
                        position=1;
                        benchmarkLong=dataList.get(k).high;
                        trades.add(1.0);

                        //先结束short position
                        Date_out.add(dataList.get(k).date);
                        Time_out.add(dataList.get(k).time);
                        Price_out.add(HH.get(k - outSampleStart));
                        Absolute_PnL.add(delta_sum);
                        Percentage_PnL.add(Absolute_PnL.get(Absolute_PnL.size() - 1)/E0);
                        delta_sum = 0.0;

                        //再开始long position
                        Date_in.add(dataList.get(k).date);
                        Time_in.add(dataList.get(k).time);
                        Price_in.add(HH.get(k - outSampleStart));
                    }
                }
                
                benchmarkShort=Math.min(dataList.get(k).low,benchmarkShort);//更新benchmarkShort为当前最低价
            }
    
            if (position== 0 && traded){} //do nothing
            
            if (position== 1 && traded){} //do nothing
            
            if (position==-1 && traded){} //do nothing
            
            //update equity
            if (k-outSampleStart != 0){
                E.add(E.get(k-outSampleStart-1)+delta);
            }
            else{
                E.set(0, E0+delta);
            }
            //calculate drawdown
            Emax=Math.max(Emax, E.get(k-outSampleStart));
            DD.add(E.get(k-outSampleStart)-Emax);
        }

        //if we still long or short futures in the last 5-min interval out of sample, we need to end our position
        if (position != 0){
            //System.out.println("Ending position");
            Date_out.add(dataList.get(outSampleEnd).date);
            Position_before_exit.add(position);
            Time_out.add(dataList.get(outSampleEnd).time);
            Price_out.add(dataList.get(outSampleEnd).close); 
            Absolute_PnL.add(delta_sum);
            Percentage_PnL.add(Absolute_PnL.get(Absolute_PnL.size() - 1)/E0);
            delta_sum = 0.0;
        }

        Double netProfit = E.get(E.size() - 1) - E0;
        System.out.println("OutSample Net Profit: " + netProfit);
        //System.out.println("InSample Max DrawDown: " + dd_max);
        outSampleResult result = new outSampleResult(Date_in, Date_out, Time_in, Time_out, Position_before_exit, 
        Price_in, Price_out, Absolute_PnL, Percentage_PnL, E, DD, recordOptParameters, startDate, endDate);
        return result;
    }

    public static List<data_i> readCSV(String csvFile) {
        String line;
        String csvSplitBy = ",";
        List<data_i> dataList = new ArrayList<>();
        boolean firstLineSkipped = false;

        try {
            BufferedReader br = new BufferedReader(new FileReader(csvFile));
            while ((line = br.readLine()) != null) {
                if (!firstLineSkipped) {
                    // Skip the first line
                    firstLineSkipped = true;
                    continue;
                }

                String[] data = line.split(csvSplitBy);
                data_i d = new data_i(data[0], data[1], Double.parseDouble(data[2]), Double.parseDouble(data[3]),Double.parseDouble(data[4]), Double.parseDouble(data[5]), Double.parseDouble(data[6]));
                dataList.add(d);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return dataList;
    }

   
    public static Double getMax(List<data_i> data){
        if (data.isEmpty()) {
            throw new IllegalArgumentException("Invalid input data or column index");
        }

        Double max = Double.MIN_VALUE;
        for (data_i d : data) {
            if (d.high > max) {
                max = d.high;
            }
        }
        return max;
    }

    public static Double getMin(List<data_i> data){
        if (data.isEmpty() ) {
            throw new IllegalArgumentException("Invalid input data or column index");
        }

        Double min = Double.MAX_VALUE;
        for (data_i d : data) {
            if (d.low < min) {
                min = d.low;
            }
        }
        return min;
    }

    public static Double getMin_double(List<Double> list) {
        if (list.isEmpty()) {
            throw new IllegalArgumentException("List is empty.");
        }

        Double min = list.get(0); // Initialize min with the first element

        // Iterate through the list to find the minimum value
        for (int i = 1; i < list.size(); i++) {
            double current = list.get(i);
            if (current < min) {
                min = current;
            }
        }

        return min;
    }

    public static int getIndex(List<data_i> list, String date){
        if (list.isEmpty()){
            throw new IllegalArgumentException("List is empty.");
        }

        int indx = Integer.MAX_VALUE;
        for (int i = 0; i < list.size(); i++){
            if (list.get(i).date.equals(date)){
                indx = i;
                return indx;
            }
        }
        return indx;
    }

    public static Double[] findOptimal2(int slpg, int PV, Double E0, List<data_i> dataList, List<Double> ChnLen, List<Double> StpPct, int inSampleStart, int inSampleEnd){
        int[] i_j = new int[2];
        Double[] optParameters = {0.0, 0.0};
        Double optRatio = -99999.0;
        int position = 0;
        Hashtable<int[], Double> inSampleResult = new Hashtable<>();

        // 1000 for channelength
        for (int i = 0; i < ChnLen.size(); i += 100){
            Double L = ChnLen.get(i);
            //System.out.println("calculating length = " + L);

            List<Double> HH = new ArrayList<>();
            List<Double> LL = new ArrayList<>();
            //get HH and LL for given length
            for (int k = inSampleStart; k <= inSampleEnd; k += 1){
                HH.add(getMax(dataList.subList(k - L.intValue(), k-1)));
                LL.add(getMin(dataList.subList(k - L.intValue(), k-1)));
            }

            //iterate through each StpPct
            for (int j = 0; j < StpPct.size(); j += 10){
                Double S = StpPct.get(j);
                //System.out.println("calculating PctStp = " + S);
                
                //set initial conditions
                position = 0;
                List<Double> E = new ArrayList<>();
                E.add(E0);
                List<Double> DD = new ArrayList<>();
                DD.add(0.0);
                List<Double> trades = new ArrayList<>();
                trades.add(0.0);
                Double Emax = E0;
                boolean traded, buy, sell, sellShort, buyLong;
                Double delta;
                Double  benchmarkLong = 0.0;
                Double benchmarkShort = 0.0;

                //running through insample and trading
                for (int k= inSampleStart; k <= inSampleEnd; k += 1){
                    traded=false;
                    delta=PV*(dataList.get(k).close-dataList.get(k-1).close)*position;
                    //System.out.println(delta);

            
                    if (position== 0){
                        //System.out.println(dataList.get(k).high);
                        //System.out.println(HH.get(k - inSampleStart - 1));
                        buy=dataList.get(k).high>=HH.get(k - inSampleStart);
                        sell=dataList.get(k).low<=LL.get(k - inSampleStart);
                        
                    
                 
                        if (buy && sell){
                            delta = -slpg+PV*(LL.get(k - inSampleStart)-HH.get(k- inSampleStart));
                            trades.add(1.0);
                        }
                        else{
                            if(buy){
                                delta = -slpg/2 + PV*(dataList.get(k).close-HH.get(k- inSampleStart));
                                position= 1;
                                traded=true;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(0.5);
                            }
                
                            if(sell){
                                //System.out.println("k: " + k);
                                //System.out.println("sell :" + sell);
                                delta = -slpg/2 - PV*(dataList.get(k).close-LL.get(k- inSampleStart));
                                //System.out.println("Delta: "+ delta);
                                position=-1;
                                traded=true;
                                benchmarkShort=dataList.get(k).low;
                                //System.out.println("benchmarkShort: " + benchmarkShort);
                                trades.add( 0.5);
                            }
                        }
                    }
                    
            
            
                    if (position== 1 && !traded){
                        sellShort=dataList.get(k).low<=LL.get(k - inSampleStart);
                        sell=dataList.get(k).low<=(benchmarkLong*(1-S));
                            
                        if(sellShort && sell){
                            //copy of sell short
                            if(sellShort){
                                    delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - inSampleStart));
                                    position=-1;
                                    benchmarkShort=dataList.get(k).low;
                                    trades.add(1.0);
                            }
                        }
                        else{
                            if(sell){
                                delta=delta-slpg/2-PV*(dataList.get(k).close-(benchmarkLong*(1-S))); //min(Open,stopPrice)
                                position=0;
                                trades.add(0.5);
                            }
                                
                            if(sellShort){
                                delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - inSampleStart)); //min(Open,LL(k))
                                position=-1;
                                benchmarkShort=dataList.get(k).low;
                                trades.add(1.0);
                            }
                        }
                            
                        benchmarkLong=Math.max(dataList.get(k).high,benchmarkLong);
                    }
            
                    if (position==-1 && !traded){
                        buyLong=dataList.get(k).high>=HH.get(k - inSampleStart);
                        buy=dataList.get(k).high>=(benchmarkShort*(1+S));
                        
                        if(buyLong && buy){
                            //copy of buyLong
                            if(buyLong){
                                delta=delta-slpg+2*PV*(dataList.get(k).close - HH.get(k - inSampleStart));
                                position=1;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(1.0);
                            }
                        }
                        else{
                            if(buy){
                                //System.out.println("k: " + k);
                                //System.out.println("short position - Buy");
                                delta=delta-slpg/2+PV*(dataList.get(k).close-(benchmarkShort*(1+S)));
                                //System.out.println("Delta: "+ delta);
                                position=0;
                                trades.add(0.5);
                            }
                            
                            if(buyLong){
                                delta=delta-slpg+2*PV*(dataList.get(k).close-HH.get(k - inSampleStart));
                                position=1;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(1.0);
                            }
                        }
                        
                        benchmarkShort=Math.min(dataList.get(k).low,benchmarkShort);
                    }
            
                    if (position== 0 && traded){} //do nothing
                    
                    if (position== 1 && traded){} //do nothing
                    
                    if (position==-1 && traded){} //do nothing
                    
                    //update equity
                    //System.out.println("Delta = " + delta);
                    //update equity
                    if (k-inSampleStart != 0){
                        E.add(E.get(k-inSampleStart-1)+delta);
                    }
                    else{
                        E.set(0, E0+delta);
                    }
                    //calculate drawdown
                    Emax=Math.max(Emax, E.get(k-inSampleStart));
                    DD.add(E.get(k-inSampleStart)-Emax);
                    //System.out.println(E.get(k - inSampleStart));
                    //System.out.println(DD.get(k - inSampleStart));
                }
                
                Double dd_max = getMin_double(DD);
                Double netProfit = E.get(E.size() - 1) - E0;
                //System.out.println("InSample Net Profit: " + netProfit);
                //System.out.println("InSample Max DrawDown: " + dd_max);
                Double netProfit_ddMax_ratio = netProfit/(-dd_max);
                //System.out.println("netProfit_ddMax_ratio: " + netProfit_ddMax_ratio);
                Double[] L_S = {L,S};
                
                // if (i_j == null) {
                //     i_j[0] = i;
                //     i_j[1] = j;  
                // }
                System.out.println("ratio: " + netProfit_ddMax_ratio);
                if (netProfit_ddMax_ratio > optRatio){
                    optRatio = netProfit_ddMax_ratio;
                    i_j[0] = i;
                    i_j[1] = j;
                    optParameters = L_S;
                    System.out.println("optimal_ratio: " + optRatio);
                }            
                //inSampleResult.put(i_j, netProfit_ddMax_ratio);

            }
        }
        System.out.println("i"  + i_j[0] );
        //现在是100 chanelength的循环
        int u = i_j[0];
        int temp = u;
        
        for (u = i_j[0]; u < temp+50; u += 10){
            Double L = ChnLen.get(u);
            //System.out.println("calculating length = " + L);

            List<Double> HH = new ArrayList<>();
            List<Double> LL = new ArrayList<>();
            //get HH and LL for given length
            for (int k = inSampleStart; k <= inSampleEnd; k += 1){
                HH.add(getMax(dataList.subList(k - L.intValue(), k-1)));
                LL.add(getMin(dataList.subList(k - L.intValue(), k-1)));
            }

            //iterate through each StpPct
            for (int j = 0; j < StpPct.size(); j += 10){
                Double S = StpPct.get(j);
                //System.out.println("calculating PctStp = " + S);
                
                //set initial conditions
                position = 0;
                List<Double> E = new ArrayList<>();
                E.add(E0);
                List<Double> DD = new ArrayList<>();
                DD.add(0.0);
                List<Double> trades = new ArrayList<>();
                trades.add(0.0);
                Double Emax = E0;
                boolean traded, buy, sell, sellShort, buyLong;
                Double delta;
                Double  benchmarkLong = 0.0;
                Double benchmarkShort = 0.0;

                //running through insample and trading
                for (int k= inSampleStart; k <= inSampleEnd; k += 1){
                    traded=false;
                    delta=PV*(dataList.get(k).close-dataList.get(k-1).close)*position;
                    //System.out.println(delta);

            
                    if (position== 0){
                        //System.out.println(dataList.get(k).high);
                        //System.out.println(HH.get(k - inSampleStart - 1));
                        buy=dataList.get(k).high>=HH.get(k - inSampleStart);
                        sell=dataList.get(k).low<=LL.get(k - inSampleStart);
                        
                    
                 
                        if (buy && sell){
                            delta = -slpg+PV*(LL.get(k - inSampleStart)-HH.get(k- inSampleStart));
                            trades.add(1.0);
                        }
                        else{
                            if(buy){
                                delta = -slpg/2 + PV*(dataList.get(k).close-HH.get(k- inSampleStart ));
                                position= 1;
                                traded=true;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(0.5);
                            }
                
                            if(sell){
                                //System.out.println("k: " + k);
                                //System.out.println("sell :" + sell);
                                delta = -slpg/2 - PV*(dataList.get(k).close-LL.get(k- inSampleStart));
                                //System.out.println("Delta: "+ delta);
                                position=-1;
                                traded=true;
                                benchmarkShort=dataList.get(k).low;
                                //System.out.println("benchmarkShort: " + benchmarkShort);
                                trades.add( 0.5);
                            }
                        }
                    }
                    
            
            
                    if (position== 1 && !traded){
                        sellShort=dataList.get(k).low<=LL.get(k - inSampleStart);
                        sell=dataList.get(k).low<=(benchmarkLong*(1-S));
                            
                        if(sellShort && sell){
                            //copy of sell short
                            if(sellShort){
                                    delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - inSampleStart));
                                    position=-1;
                                    benchmarkShort=dataList.get(k).low;
                                    trades.add(1.0);
                            }
                        }
                        else{
                            if(sell){
                                delta=delta-slpg/2-PV*(dataList.get(k).close-(benchmarkLong*(1-S))); //min(Open,stopPrice)
                                position=0;
                                trades.add(0.5);
                            }
                                
                            if(sellShort){
                                delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - inSampleStart)); //min(Open,LL(k))
                                position=-1;
                                benchmarkShort=dataList.get(k).low;
                                trades.add(1.0);
                            }
                        }
                            
                        benchmarkLong=Math.max(dataList.get(k).high,benchmarkLong);
                    }
            
                    if (position==-1 && !traded){
                        buyLong=dataList.get(k).high>=HH.get(k - inSampleStart);
                        buy=dataList.get(k).high>=(benchmarkShort*(1+S));
                        
                        if(buyLong && buy){
                            //copy of buyLong
                            if(buyLong){
                                delta=delta-slpg+2*PV*(dataList.get(k).close - HH.get(k - inSampleStart));
                                position=1;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(1.0);
                            }
                        }
                        else{
                            if(buy){
                                //System.out.println("k: " + k);
                                //System.out.println("short position - Buy");
                                delta=delta-slpg/2+PV*(dataList.get(k).close-(benchmarkShort*(1+S)));
                                //System.out.println("Delta: "+ delta);
                                position=0;
                                trades.add(0.5);
                            }
                            
                            if(buyLong){
                                delta=delta-slpg+2*PV*(dataList.get(k).close-HH.get(k - inSampleStart));
                                position=1;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(1.0);
                            }
                        }
                        
                        benchmarkShort=Math.min(dataList.get(k).low,benchmarkShort);
                    }
            
                    if (position== 0 && traded){} //do nothing
                    
                    if (position== 1 && traded){} //do nothing
                    
                    if (position==-1 && traded){} //do nothing
                    
                    //update equity
                    //System.out.println("Delta = " + delta);
                    //update equity
                    if (k-inSampleStart != 0){
                        E.add(E.get(k-inSampleStart-1)+delta);
                    }
                    else{
                        E.set(0, E0+delta);
                    } 
                    //calculate drawdown
                    Emax=Math.max(Emax, E.get(k-inSampleStart));
                    DD.add(E.get(k-inSampleStart)-Emax);
                    //System.out.println(E.get(k - inSampleStart));
                    //System.out.println(DD.get(k - inSampleStart));
                }
                
                Double dd_max = getMin_double(DD);
                Double netProfit = E.get(E.size() - 1) - E0;
                //System.out.println("InSample Net Profit: " + netProfit);
                //System.out.println("InSample Max DrawDown: " + dd_max);
                Double netProfit_ddMax_ratio = netProfit/(-dd_max);
                //System.out.println("netProfit_ddMax_ratio: " + netProfit_ddMax_ratio);
                Double[] L_S = {L,S};
                
                //int[] u_j = {u, j};
                if (netProfit_ddMax_ratio > optRatio){
                    optRatio = netProfit_ddMax_ratio;
                    optParameters = L_S;
                    i_j[0] = u;
                    i_j[1] = j;
                }            
                //inSampleResult.put(u_j, netProfit_ddMax_ratio);

            }
        }
        //现在是10 chanelength的循环
        int y = i_j[0];
        temp = y;
        for (y = i_j[0]; y < temp+5; y += 1){
            Double L = ChnLen.get(y);
            System.out.println("calculating length = " + L);

            List<Double> HH = new ArrayList<>();
            List<Double> LL = new ArrayList<>();
            //get HH and LL for given length
            for (int k = inSampleStart; k <= inSampleEnd; k += 1){
                HH.add(getMax(dataList.subList(k - L.intValue(), k-1)));
                LL.add(getMin(dataList.subList(k - L.intValue(), k-1)));
            }

            //iterate through each StpPct
            for (int j = 0; j < StpPct.size(); j += 1){
                Double S = StpPct.get(j);
                //System.out.println("calculating PctStp = " + S);
                
                //set initial conditions
                position = 0;
                List<Double> E = new ArrayList<>();
                E.add(E0);
                List<Double> DD = new ArrayList<>();
                DD.add(0.0);
                List<Double> trades = new ArrayList<>();
                trades.add(0.0);
                Double Emax = E0;
                boolean traded, buy, sell, sellShort, buyLong;
                Double delta;
                Double  benchmarkLong = 0.0;
                Double benchmarkShort = 0.0;

                //running through insample and trading
                for (int k= inSampleStart; k <= inSampleEnd; k += 1){
                    traded=false;
                    delta=PV*(dataList.get(k).close-dataList.get(k-1).close)*position;
                    //System.out.println(delta);

            
                    if (position== 0){
                        //System.out.println(dataList.get(k).high);
                        //System.out.println(HH.get(k - inSampleStart - 1));
                        buy=dataList.get(k).high>=HH.get(k - inSampleStart);
                        sell=dataList.get(k).low<=LL.get(k - inSampleStart);
                        
                    
                 
                        if (buy && sell){
                            delta = -slpg+PV*(LL.get(k - inSampleStart)-HH.get(k- inSampleStart));
                            trades.add(1.0);
                        }
                        else{
                            if(buy){
                                delta = -slpg/2 + PV*(dataList.get(k).close-HH.get(k- inSampleStart ));
                                position= 1;
                                traded=true;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(0.5);
                            }
                
                            if(sell){
                                //System.out.println("k: " + k);
                                //System.out.println("sell :" + sell);
                                delta = -slpg/2 - PV*(dataList.get(k).close-LL.get(k- inSampleStart));
                                //System.out.println("Delta: "+ delta);
                                position=-1;
                                traded=true;
                                benchmarkShort=dataList.get(k).low;
                                //System.out.println("benchmarkShort: " + benchmarkShort);
                                trades.add( 0.5);
                            }
                        }
                    }
                    
            
            
                    if (position== 1 && !traded){
                        sellShort=dataList.get(k).low<=LL.get(k - inSampleStart);
                        sell=dataList.get(k).low<=(benchmarkLong*(1-S));
                            
                        if(sellShort && sell){
                            //copy of sell short
                            if(sellShort){
                                    delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - inSampleStart));
                                    position=-1;
                                    benchmarkShort=dataList.get(k).low;
                                    trades.add(1.0);
                            }
                        }
                        else{
                            if(sell){
                                delta=delta-slpg/2-PV*(dataList.get(k).close-(benchmarkLong*(1-S))); //min(Open,stopPrice)
                                position=0;
                                trades.add(0.5);
                            }
                                
                            if(sellShort){
                                delta=delta-slpg-2*PV*(dataList.get(k).close-LL.get(k - inSampleStart)); //min(Open,LL(k))
                                position=-1;
                                benchmarkShort=dataList.get(k).low;
                                trades.add(1.0);
                            }
                        }
                            
                        benchmarkLong=Math.max(dataList.get(k).high,benchmarkLong);
                    }
            
                    if (position==-1 && !traded){
                        buyLong=dataList.get(k).high>=HH.get(k - inSampleStart);
                        buy=dataList.get(k).high>=(benchmarkShort*(1+S));
                        
                        if(buyLong && buy){
                            //copy of buyLong
                            if(buyLong){
                                delta=delta-slpg+2*PV*(dataList.get(k).close - HH.get(k - inSampleStart));
                                position=1;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(1.0);
                            }
                        }
                        else{
                            if(buy){
                                //System.out.println("k: " + k);
                                //System.out.println("short position - Buy");
                                delta=delta-slpg/2+PV*(dataList.get(k).close-(benchmarkShort*(1+S)));
                                //System.out.println("Delta: "+ delta);
                                position=0;
                                trades.add(0.5);
                            }
                            
                            if(buyLong){
                                delta=delta-slpg+2*PV*(dataList.get(k).close-HH.get(k - inSampleStart));
                                position=1;
                                benchmarkLong=dataList.get(k).high;
                                trades.add(1.0);
                            }
                        }
                        
                        benchmarkShort=Math.min(dataList.get(k).low,benchmarkShort);
                    }
            
                    if (position== 0 && traded){} //do nothing
                    
                    if (position== 1 && traded){} //do nothing
                    
                    if (position==-1 && traded){} //do nothing
                    
                    //update equity
                    //System.out.println("Delta = " + delta);
                    //update equity
                    if (k-inSampleStart != 0){
                        E.add(E.get(k-inSampleStart-1)+delta);
                    }
                    else{
                        E.set(0, E0+delta);
                    } 
                    //calculate drawdown
                    Emax=Math.max(Emax, E.get(k-inSampleStart));
                    DD.add(E.get(k-inSampleStart)-Emax);
                    //System.out.println(E.get(k - inSampleStart));
                    //System.out.println(DD.get(k - inSampleStart));
                }
                
                Double dd_max = getMin_double(DD);
                Double netProfit = E.get(E.size() - 1) - E0;
                //System.out.println("InSample Net Profit: " + netProfit);
                //System.out.println("InSample Max DrawDown: " + dd_max);
                Double netProfit_ddMax_ratio = netProfit/(-dd_max);
                //System.out.println("netProfit_ddMax_ratio: " + netProfit_ddMax_ratio);
                Double[] L_S = {L,S};
            
                //int[] y_j = {y, j};
                if (netProfit_ddMax_ratio > optRatio){
                    optRatio = netProfit_ddMax_ratio;
                    optParameters = L_S;
                    System.out.println("optimal L_S" +": "+ L +','+ S);
                    i_j[0] = y;
                    i_j[1] = j;
                }            
                inSampleResult.put(i_j, netProfit_ddMax_ratio);
        
            }
        
        }
        return optParameters;
    }
}



