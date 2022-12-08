/* The following code will model a trading stock market simulator will simulate the performance of certain portfolios.
** The project is divided into 2 simulations: the market simulation and the investor simulation. Interactions between
** both will later be included, creating cause and effect relationships.
** 
** Market Simulation: Creation of stochastic process through mathematical model, in which performance of each share
** is represented by a pseudo-random function that evolves over time. Trend of shares are characterised by statistical
** parameters that will determine their volatility over time. [Previously know data set, random function for unknown
** data points after t=0].
**
** Investor Simulation: Involves management of stock wallets, with each able to follow invest strategies. First stage
** [current], market shares will be assigned at beginning of simulation and stay assigned until end of simulation.
** Later, wallets will be able to buy and sell market shares during the simulation by following the strategies.
** 
** Data Set: Using S&P500 datapoints from 02022013-05022018 [5 years], we obtained all stocks available, calculating
** their statistical parameters.
*/

#include <iostream> //
#include <sstream> //
#include <fstream> //
#include <vector> //
#include <string> //
#include <unordered_map> //
#include <cmath> //
#include <memory>
#include <algorithm> //
#include <array> 
#include <random> //
#include <functional> 

// MATHEMATICAL FUNCTIONS
void MultiplyVectorByScalar(std::vector<int> &v, int k){
	for_each(v.begin(), v.end(), [k](int &c){ c *= k; });
}

double scalar_prod(std::vector<double> v1, std::vector<double> v2){
    
    double s_prod = 0.0;

    for(int i=1;i < v1.size(); ++i){
        s_prod = (v1[i]*v2[i]) + s_prod;
    }
    return s_prod;
}



//***************************************************************************************************************
//************************************************USER INPUT*****************************************************
//***************************************************************************************************************

// Function that allows the user to input funds to invest in the market
double inputFunds() 
{
    double funds;
    std::cout << "Enter the funds to invest in the market: ";
    std::cin >> funds;
    return funds;
}

// Define a function that allows the user to choose a strategy from a list of predefined strategies
int chooseStrategy(std::vector<std::string> strategies)
{
    std::cout << "Available Asset Allocation Strategies: " << std::endl;      // Print the available strategies
    for (int i = 0; i < strategies.size(); ++i)
    {
        std::cout << i + 1 << ". " << strategies[i] << std::endl;
    }

    std::cout << "Choose an Asset Allocation Strategy: "; // Prompt the user to choose a strategy
    int index;                          // Read the user's input
    std::cin >> index;

    if (index < 1 || index > strategies.size()) // Check user's input validity
    {
        std::cout << "\n \n";
        std::cout << "Please repeat ID selection" << std::endl; // Error message
        return chooseStrategy(strategies);  // Recursively call chooseStrategy function to get a valid strategy
    }
    return index;       // Return the chosen strategy
}

std::vector<std::string> inputCompanies(){

    std::string share_names;
    std::cout << "\n Enter the names of the shares you want to use, without spaces and separated by a comma: ";
    std::cin >> share_names;

    // Split the input string on the comma separator
    std::vector<std::string> selected_shares;
    std::string current_share;
    for (char c : share_names) {
        if (c == ',') {
        selected_shares.push_back(current_share);
        current_share.clear();
        } else {
        current_share += c;
        }
    }

    if (!current_share.empty()) {                       // Initialize vector of strings with no current shares
        selected_shares.push_back(current_share);
    }

    return selected_shares;
}


int inputTime() 
{
    std::cout << "How long do you want the simulation to go for? (days) ";
    int time_f;                                             // Read the user's input
    std::cin >> time_f;
    return time_f;
}


// A class representing an asset
class Assets {
public:
  std::string name; //Asset Name. To be later completed with more parameters of the class
  Assets(std::string name) : name(name) {}  // Assets class constructor
};

// A class representing a share, inherited from Assets
class Shares : public Assets {
public:
  double start_value;       //Intial Market Value of Share
  double variance;          //Statistic Parameter 1
  double price_change;      //Statistic Parameter 2
  double expected_return;   //Statistic Parameter 3

  Shares(std::string name, double start_value, double variance, 
  double price_change, double expected_return) : Assets(name), 
  start_value(start_value), variance(variance), price_change(price_change), 
  expected_return(expected_return) {}                                           // Shares class constructor
};




//***************************************************************************************************************
//*********************************************MARKET SIMULATION*************************************************
//***************************************************************************************************************



// Structure with the S&P500 data
struct SP500Data
{
    std::vector<std::string> company;

    // The following variables change over the date, so we should define a vector of vectors
    std::vector<std::vector<double>> close;

    //At the moment, our expected return & share price stimation only uses dayly close price  
    // Uncoment if nedded for new fomulation   
    // std::vector<std::vector<double>> volume;
    // std::vector<std::vector<double>> open;
    // std::vector<std::vector<double>> high;
    // std::vector<std::vector<double>> low;
    // std::vector<std::vector<double>> volume;
};




std::unique_ptr<SP500Data> read_from_csv(std::string filename="all_stocks_5yr.csv")
{
    // Read the data from the csv 
    std::ifstream database(filename);

    std::unique_ptr<SP500Data> data = std::make_unique<SP500Data>();

    std::string date;
    double close;

    // Even this variables are not saved in the structure, we should decare them to avoid read-location errors 
    // csv strings differ in char numbers, so at the moment this is the easiest way to load data.
    //TODO: Optimice this 
    double open;
    double high;
    double low;
    double volume;

    std::string company;
    std::string curr_company = "";
    auto company_id = -1;

    std::string line; 

    std::getline(database, line);  // Header line remove


    while (std::getline(database, line)){

            std::replace(line.begin(), line.end(), ',', ' ');

            // stringstream breaks on spaces, check variable declarations to avoid errors (string, int, float ...)
            std::stringstream line_stream(line);
            

            line_stream >> date >> open >> high >> low >> close >> volume >> company;
            // csv header names: date,open,high,low,close,volume,Name

            if( curr_company != company ){

                curr_company = company;
                data->company.push_back(company);
                data->close.push_back(std::vector<double>{});
            //     data->volume.push_back(std::vector<double>{});
            //     data->open.push_back(std::vector<double>{});
            //     data->high.push_back(std::vector<double>{});
            //     data->low.push_back(std::vector<double>{});
                ++company_id;
            }
            
            data->close[company_id].push_back(close);
    }
    return data;
}

//Population of random parameters [min,max] dimmension 'num' double
void random_vec( std::vector<double>& num) {

    auto seed = 0;
    auto min = -1.0;
    auto max = 1.0;

    std::mt19937 gen(seed);

    // uniform_int_distribution: use the mersenne twister
    // engine to generate a uniform random distribution over (min, max)
    std::uniform_real_distribution<double> unif_distrib(min, max);
   // 3. solution
   for( auto& elem : num){
       elem = unif_distrib(gen);
   }
}

/*Compute gbm estimation over a time 't'. 
Imput: Share parameters double[3]
    S    = Stock price
    μ    = Expected return
    σ    = Standard deviation of returns
    Δt   = Elapsed time period
Output: S_t= Change in stock price at time t */
std::vector<double> gbm_2( double S, double nu, double dev, int t_fin) {

    double gbm_operator;
    double delta_S;
    
    int delta_t = 1; // Time step for iterations. ex = one day. Must be INT

    std::vector<double> S_trend;

    std::vector<double> eps(t_fin,0);
    random_vec(eps);
    
    for(int t = 0; t <= t_fin; t++)
    {
        gbm_operator = (nu - dev*dev/2)*delta_t + dev*eps[t]*sqrt(delta_t);
        
        S = S * exp(gbm_operator); // Compute S in t+1

        S_trend.push_back(S); // Saves the stock price change
    }

    return S_trend;
}

double standar_dev(std::vector<double> data){

    double mean;
    double opr=0.0;
    double standar_dev;

    // Higher share price Iterator for normalice the standar deviation
    std::vector<double>::iterator max_value;
    // Data dimension (double)
    double dim = double(data.size());

    mean = std::accumulate( data.begin(), data.end(), 0 ) / dim;

    for(int i = 0; i<data.size(); ++i){
       opr = (data[i]-mean)*(data[i]-mean) + opr;
    };

    
    max_value = std::max_element(data.begin(), data.end());

    standar_dev = sqrt(opr / dim) / *max_value ;

    return standar_dev;
};

double expected_return(std::vector<double> data){

    // historic return average  
    std::vector<double> return_prcn;
    double expected_return;
    double dim = double(data.size());


    for(int i = 0; i<data.size()-1; ++i){
       return_prcn.push_back( (data[i+1] - data[i]) / data[i] );
    };

    expected_return = std::accumulate( return_prcn.begin(), return_prcn.end(), 0.0 ) / dim;
    
    return expected_return;
};


// Search the names of the selected companies in bigger string vectors (from structure) and return the 
// company ID of each name
void name_to_id( std::vector<int>& v_id, std::vector<std::string>& names, std::vector<std::string>& nameslist)
{   
    // std::unique_ptr<SP500Data> data_frame;
    // // data_frame->company
    // std::vector<std::string> nameslist = data_frame->company;

    names = inputCompanies();// Load the name of the companies

    for(auto& elem : names){ // Short list loop 
        for( int i = 0; i < nameslist.size(); ++i ) // Long list loop
        {
            if( elem == nameslist[i]){
                v_id.push_back(i); // Compute the companies ID vector in order
            }
        }
    }

    if (v_id.size() < 1 ||  v_id.size() != names.size() ) // Check user's input validity
    {
        std::cout << "\n \n";
        std::cout << "Please repeat Companies selection" << std::endl; // Error message
        
        v_id.clear();// Clear the vector values
        names.clear();
        return name_to_id(v_id, names, nameslist);  // Recursively call chooseStrategy function to get a valid strategy
    }

   
}



// Compute the market simulation and print the results for each selected share
void shares_sim(std::vector<int> selected_company_ids, std::vector<std::vector<double>> close,
                std::vector<std::string> selected_shares, int t_f, std::vector<double> vec_fps){
                    
    double expected_r;
    double dev; 
    
    double portfolio_returns=0;

    for (int i = 0; i < selected_company_ids.size(); i++){
        std::cout << "************************************************** "  << std::endl ;
        std::cout << "Company Name: " << selected_shares[i]
                            << "    Company ID: " << selected_company_ids[i] << std::endl ;

        expected_r = expected_return(close[selected_company_ids[i]]);
        dev = standar_dev(close[selected_company_ids[i]]);
        std::cout << "expected_return: " << expected_r << "    standar_dev: " << dev << std::endl ;

        std::vector<double> S_trend;
        
        S_trend = gbm_2(close[selected_company_ids[i]].back(), expected_r, dev, t_f);
        

        for (int i = 0; i < S_trend.size(); i++)
        {
            std::cout << " step " << i+1 << ": "<< S_trend[i] << std::endl ;
        }
        
        
        std::cout << "Initial price: " << close[selected_company_ids[i]].back()
                     << " --> Final price: " << S_trend.back() << std::endl ;
        
        double partial_return = vec_fps[i] * (S_trend.back()-close[selected_company_ids[i]].back()) / close[selected_company_ids[i]].back(); 
        std::cout << "These are the partial returns: " << partial_return << std::endl;
        portfolio_returns += partial_return;
    }
    std::cout << "Portfolio Returns :\t" << portfolio_returns << std::endl;
}


// Search the stock company names into the preloaded defined groups 
std::vector<int> name_to_category(std::vector<std::string> names)
{
    // Company names defined directly at the function to avoid function load and output (Researh & confirm runtime best option)
    std::vector<std::string> low = {"MA", "MAR", "NOV", "A", "AAP", "ABBV", "ABC", "ABT", "ACN", "ADM", "ADP", "ADS", "AEE", "AEP", "AES", "AET", "AFL", "AGN", "AIG", "AIV", "AIZ", "AJG", "ALB", "ALK", "ALL", "ALLE", "ALXN", "AMAT", "AMD", "AME", "AMG", "AMP", "AMT", "AMZN", "ANDV", "ANSS", "ANTM", "AON", "APA", "APC", "APD", "APH", "APTV", "ARE", "ARNC", "ATVI", "AVB", "AVY", "AWK", "AXP", "AYI", "AZO", "BA", "BAC", "BAX", "BBT", "BBY", "BDX", "BF.B", "BHGE", "BK", "BLK", "BLL", "BMY", "BRK.B", "BSX", "BWA", "BXP", "C", "CA", "CAG", "CAH", "CAT", "CB", "CBG", "CBOE", "CBS", "CCI", "CCL", "CDNS", "CERN", "CF", "CFG", "CHD", "CHK", "CHRW", "CI", "CINF", "CL", "CLX", "CMA", "CME", "CMG", "CMI", "CMS", "CNC", "CNP", "COF", "COL", "COO", "COP", "COTY", "CPB", "CSRA", "CSX", "CTAS", "CTL", "CTXS", "CVS", "CVX", "CXO", "D", "DAL", "DE", "DFS", "DG", "DGX", "DHI", "DHR", "DIS", "DISCK", "DISH", "DLR", "DLTR", "DOV", "DPS", "DRE", "DRI", "DTE", "DUK", "DVA", "DVN", "DWDP", "DXC", "ECL", "ED", "EFX", "EIX", "EL", "EMN", "EMR", "EOG", "EQIX", "EQR", "EQT", "ES", "ESS", "ETFC", "ETN", "ETR", "EVHC", "EW", "EXC", "EXPD", "EXPE", "EXR", "F", "FAST", "FBHS", "FCX", "FDX", "FE", "FIS", "FISV", "FITB", "FL", "FLIR", "FLR", "FMC", "FOX", "FRT", "FTI", "FTV", "GD", "GE", "GGP", "GIS", "GLW", "GM", "GPC", "GPN", "GPS", "GS", "GWW", "HAL", "HAS", "HCA", "HCN", "HCP", "HD", "HES", "HIG", "HII", "HLT", "HOG", "HOLX", "HON", "HP", "HPE", "HPQ", "HRB", "HRL", "HRS", "HST", "HSY", "HUM", "IBM", "ICE", "IFF", "ILMN", "INFO", "INTU", "IP", "IPG", "IQV", "IR", "IRM", "ISRG", "IT", "ITW", "IVZ", "JCI", "JEC", "JNJ", "JNPR", "JPM", "JWN", "K", "KEY", "KHC", "KIM", "KLAC", "KMB", "KMI", "KMX", "KO", "KORS", "KR", "KSS", "KSU", "LB", "LEG", "LEN", "LH", "LKQ", "LLL", "LLY", "LMT", "LNC", "LNT", "LOW", "LRCX", "LUK", "LYB", "M", "MAA", "MAC", "MAS", "MAT", "MCD", "MCHP", "MCK", "MCO", "MDT", "MET", "MGM", "MHK", "MKC", "MLM", "MMC", "MMM", "MNST", "MO", "MON", "MOS", "MPC", "MRK", "MRO", "MS", "MSI", "MTB", "MTD", "NAVI", "NCLH", "NDAQ", "NEE", "NEM", "NFLX", "NFX", "NI", "NKE", "NLSN", "NOC", "NRG", "NSC", "NTAP", "NTRS", "NUE", "NVDA", "NWL", "NWS", "OKE", "OMC", "ORCL", "ORLY", "OXY", "PAYX", "PCAR", "PCLN", "PDCO", "PEG", "PEP", "PFE", "PFG", "PG", "PGR", "PH", "PHM", "PKG", "PKI", "PLD", "PM", "PNC", "PNR", "PNW", "PPG", "PRGO", "PRU", "PSX", "PVH", "PWR", "PX", "PXD", "PYPL", "QRVO", "RCL", "RE", "REG", "RF", "RHI", "RHT", "RJF", "RL", "RMD", "ROK", "ROP", "RRC", "RSG", "RTN", "SCG", "SCHW", "SEE", "SHW", "SIG", "SJM", "SLB", "SLG", "SNA", "SNI", "SO", "SPG", "SPGI", "SRCL", "SRE", "STI", "STT", "STZ", "SWK", "SWKS", "SYF", "SYK", "SYY", "TAP", "TDG", "TEL", "TGT", "TIF", "TJX", "TMO", "TPR", "TRIP", "TROW", "TRV", "TSCO", "TSN", "TSS", "TXN", "TXT", "UA", "UAL", "UDR", "UHS", "ULTA", "UNH", "UNM", "UNP", "UPS", "URI", "USB", "UTX", "VAR", "VIAB", "VLO", "VMC", "VNO", "VRSK", "VRSN", "VRTX", "VTR", "VZ", "WAT", "WBA", "WEC", "WFC", "WHR", "WLTW", "WM", "WMB", "WMT", "WRK", "WU", "WY", "WYN", "WYNN", "XEC"};
    std::vector<std::string> mid = {"AAPL", "ADBE", "ADI", "AKAM", "ALGN", "AMGN", "AOS", "AVGO", "CELG", "CHTR", "CMCSA", "COG", "COST", "CTSH", "DISCA", "EA", "EBAY", "ESRX", "FFIV", "FOXA", "GOOG", "GOOGL", "GT", "HBAN", "HBI", "HSIC", "IDXX", "INCY", "JBHT", "MSFT", "MYL", "NBL", "PBCT", "QCOM", "REGN", "SNPS", "STX", "SYMC", "T", "TMK", "TWX", "UAA", "V", "VFC"};
    std::vector<std::string> high = {"AAL", "ADSK", "BEN", "BHF", "BIIB", "CRM", "CSCO", "FB", "FLS", "GILD", "GRMN", "INTC", "L", "LUV", "MDLZ", "MU", "NWSA", "O", "PCG", "PPL", "PSA", "ROST", "SBAC", "SBUX", "WDC"};

    int category;

    std::vector<int> v_category;
    for(auto& elem : names){
        for( int i = 0; i < low.size(); ++i  )
        {
            if( elem == low[i]){
                category = 1;
                v_category.push_back(category);
            }
        }
        for( int i = 0; i < mid.size(); ++i  )
        {
            if( elem == mid[i]){
                category = 2;
                v_category.push_back(category);
            }
        }
        for( int i = 0; i < high.size(); ++i  )
        {
            if( elem == high[i]){
                category = 3;
                v_category.push_back(category);
            }
        }
    }
    return v_category;
}

//Counts the number of shares in a category
std::vector<double> funds_per_share(std::vector<int> vec_strat, double total_funds, int index){

    double opr;
    std::vector<std::vector<double>> strategy = {{0.6, 0.3, 0.1}, {0.2, 0.6, 0.2},{0.1, 0.1, 0.8}};
    std::vector<double> selected_strategy;
    std::vector<double> vec_fps(vec_strat.begin(), vec_strat.end());
    std::vector<int> vec_num_items;

    selected_strategy = strategy[index-1]; 
    
    
    //determine how many integers match a target value.
    for (const int target: {1, 2, 3})
    {
        const long int num_items = std::count(vec_strat.cbegin(), vec_strat.cend(), target);
        std::cout << "number: " << target << ", count: " << num_items << '\n';

        vec_num_items.push_back(num_items);
        opr = selected_strategy[target-1]/double(num_items) * total_funds;
        std::replace(vec_fps.begin(), vec_fps.end(), double(target), opr); 
    }

    // // Check if any of the values in the vector are 0
    // std::unique_ptr<SP500Data> data_frame;

    for (int value : vec_num_items) {
        if (value == 0) {
        // Ask the user to input new values if a 0 is found
        std::cout << "You haven't followed the strategy, meaning that you didn't choose a minimum of 1 or more stocks per category. \n Reinitiate to run properly, or continue to obtain a different result. \n" << std::endl;//name_to_id(data_frame->company,std::vector<int> selected_company_ids, selected_shares); // Call function
        break;
        }
    }


    return vec_fps;
}



int main()
{
    
    // Pointer declaration
    std::unique_ptr<SP500Data> data_frame;

    data_frame = read_from_csv();

    char rerun;
    do{
        std::string text = "*************************************************************************************\n"
                        "************************* WELCOME TO THE MARKET SIMULATOR ***************************\n"
                        "*************************************************************************************\n"
                        "** Find this tool as a market simulator that utilizes real world data of the SP500 **\n"
                        "* company stocks to create an investment performance analysis. Different strategies *\n"
                        "* are defined, from which you will be able to choose what assets and funds you want *\n"
                        "* to analyse. ***********************************************************************\n"
                        "*************************************************************************************\n"
                        "* Authors: Guerrero Hernandez, Anibal || Zuniga Martinez, Ignacio *******************\n"
                        "*************************************************************************************\n"
                        "* Latest Version: 08/12/2022 ********************************************************\n"
                        "*************************************************************************************\n"
                        "*************************************************************************************\n";
        std::cout << text << std::endl;
        
        double total_funds = inputFunds();                                     // Call the inputFunds function to get the user's input    
        std::cout << "Funds to invest: " << total_funds << "$ \n"<< std::endl; // Print the input funds

        std::vector<std::string> strategies = {"[Low Risk Strategy] \t Playing It Safe:\n 60\% Low Risk Assets \t|| 30\% Mid Risk Assets\t|| 10\% High Risk Assets\n", 
        "[Mid Risk Strategy] \t In for the Long Run:\n 20\% Low Risk Assets \t|| 60\% Mid Risk Assets\t|| 20\% High Risk Assets\n", 
        "[High Risk Strategy] \t Risk Off:\n 10\% Low Risk Assets \t|| 10\% Mid Risk Assets\t|| 80\% High Risk Assets\n"};     // Define a list of available strategies
        int strategy = chooseStrategy(strategies);                                                      // Call the chooseStrategy function to get the user's choice
        std::cout << "Chosen Strategy with ID = " << strategy << std::endl;    // Print the chosen strategy

        // Print the list of shares for each risk category
        std::cout << "Would you like to see all the stocks from the SP500 & their corresponding risk categories? (Y/N)" << std::endl;

        char response_see_stocks;
        std::cin >> response_see_stocks;
        std::string low_risk_shares = "MA, MAR, NOV, A, AAP, ABBV, ABC, ABT, ACN, ADM, ADP, ADS, AEE, AEP, AES, AET, AFL, AGN, AIG, AIV, AIZ, AJG, ALB, ALK, ALL, ALLE, ALXN, AMAT, AMD, AME, AMG, AMP, AMT, AMZN, ANDV, ANSS, ANTM, AON, APA, APC, APD, APH, APTV, ARE, ARNC, ATVI, AVB, AVY, AWK, AXP, AYI, AZO, BA, BAC, BAX, BBT, BBY, BDX, BF.B, BHGE, BK, BLK, BLL, BMY, BRK.B, BSX, BWA, BXP, C, CA, CAG, CAH, CAT, CB, CBG, CBOE, CBS, CCI, CCL, CDNS, CERN, CF, CFG, CHD, CHK, CHRW, CI, CINF, CL, CLX, CMA, CME, CMG, CMI, CMS, CNC, CNP, COF, COL, COO, COP, COTY, CPB, CSRA, CSX, CTAS, CTL, CTXS, CVS, CVX, CXO, D, DAL, DE, DFS, DG, DGX, DHI, DHR, DIS, DISCK, DISH, DLR, DLTR, DOV, DPS, DRE, DRI, DTE, DUK, DVA, DVN, DWDP, DXC, ECL, ED, EFX, EIX, EL, EMN, EMR, EOG, EQIX, EQR, EQT, ES, ESS, ETFC, ETN, ETR, EVHC, EW, EXC, EXPD, EXPE, EXR, F, FAST, FBHS, FCX, FDX, FE, FIS, FISV, FITB, FL, FLIR, FLR, FMC, FOX, FRT, FTI, FTV, GD, GE, GGP, GIS, GLW, GM, GPC, GPN, GPS, GS, GWW, HAL, HAS, HCA, HCN, HCP, HD, HES, HIG, HII, HLT, HOG, HOLX, HON, HP, HPE, HPQ, HRB, HRL, HRS, HST, HSY, HUM, IBM, ICE, IFF, ILMN, INFO, INTU, IP, IPG, IQV, IR, IRM, ISRG, IT, ITW, IVZ, JCI, JEC, JNJ, JNPR, JPM, JWN, K, KEY, KHC, KIM, KLAC, KMB, KMI, KMX, KO, KORS, KR, KSS, KSU, LB, LEG, LEN, LH, LKQ, LLL, LLY, LMT, LNC, LNT, LOW, LRCX, LUK, LYB, M, MAA, MAC, MAS, MAT, MCD, MCHP, MCK, MCO, MDT, MET, MGM, MHK, MKC, MLM, MMC, MMM, MNST, MO, MON, MOS, MPC, MRK, MRO, MS, MSI, MTB, MTD, NAVI, NCLH, NDAQ, NEE, NEM, NFLX, NFX, NI, NKE, NLSN, NOC, NRG, NSC, NTAP, NTRS, NUE, NVDA, NWL, NWS, OKE, OMC, ORCL, ORLY, OXY, PAYX, PCAR, PCLN, PDCO, PEG, PEP, PFE, PFG, PG, PGR, PH, PHM, PKG, PKI, PLD, PM, PNC, PNR, PNW, PPG, PRGO, PRU, PSX, PVH, PWR, PX, PXD, PYPL, QRVO, RCL, RE, REG, RF, RHI, RHT, RJF, RL, RMD, ROK, ROP, RRC, RSG, RTN, SCG, SCHW, SEE, SHW, SIG, SJM, SLB, SLG, SNA, SNI, SO, SPG, SPGI, SRCL, SRE, STI, STT, STZ, SWK, SWKS, SYF, SYK, SYY, TAP, TDG, TEL, TGT, TIF, TJX, TMO, TPR, TRIP, TROW, TRV, TSCO, TSN, TSS, TXN, TXT, UA, UAL, UDR, UHS, ULTA, UNH, UNM, UNP, UPS, URI, USB, UTX, VAR, VIAB, VLO, VMC, VNO, VRSK, VRSN, VRTX, VTR, VZ, WAT, WBA, WEC, WFC, WHR, WLTW, WM, WMB, WMT, WRK, WU, WY, WYN, WYNN, XEC";
        std::string mid_risk_shares = "AAPL, ADBE, ADI, AKAM, ALGN, AMGN, AOS, AVGO, CELG, CHTR, CMCSA, COG, COST, CTSH, DISCA, EA, EBAY, ESRX, FFIV, FOXA, GOOG, GOOGL, GT, HBAN, HBI, HSIC, IDXX, INCY, JBHT, MSFT, MYL, NBL, PBCT, QCOM, REGN, SNPS, STX, SYMC, T, TMK, TWX, UAA, V, VFC";
        std::string high_risk_shares = "AAL, ADSK, BEN, BHF, BIIB, CRM, CSCO, FB, FLS, GILD, GRMN, INTC, L, LUV, MDLZ, MU, NWSA, O, PCG, PPL, PSA, ROST, SBAC, SBUX, WDC";

        if (response_see_stocks == 'Y'){
            std::cout << "Available Assets to choose from: \n" << std::endl;
            std::cout << "Low Risk Shares: \n" << low_risk_shares << std::endl;    
            std::cout << "Mid Risk Shares: \n" << mid_risk_shares << std::endl; 
            std::cout << "High Risk Shares: \n" << high_risk_shares << std::endl;
        }
        else{}

        
        std::vector<int> selected_company_ids;
        std::vector<std::string> selected_shares;
        system("pause");

        // Ask for the companys names and output the a vector with the diferents ID
        name_to_id(selected_company_ids, selected_shares, data_frame->company); 

        std::vector<int> s_com_category = name_to_category(selected_shares);
        std::vector<double> vec_fps = funds_per_share(s_com_category, total_funds, strategy);

        for (int i = 0; i < selected_shares.size(); i++)
        {
            std::cout << "selected_shares " << selected_shares[i] <<std::endl;
            std::cout << "selected_company_ids "<< selected_company_ids[i] << std::endl;
            std::cout << "fps "<< vec_fps[i] << std::endl;
        }



        
        system("pause");
        
        // Ask for the time that the simulation would run
        int t_f;
        t_f = inputTime(); 
        
        


        
        for (int i = 0; i < selected_shares.size(); i++)
        {
            std::cout << "selected_shares " << selected_shares[i] <<std::endl;
            std::cout << "selected_company_ids "<< selected_company_ids[i] << std::endl;
            std::cout << "fps "<< vec_fps[i] << std::endl;
        }

        // Run the market simulation and print the values
        std::cout << "With the current assets and strategy, your portfolio is: " << std::endl;
        shares_sim(selected_company_ids, data_frame->close, selected_shares, t_f, vec_fps);

        
        std::cout << "DONE! Would you want to run the program again? Write Y/N\n" << std::endl;
        std::cin >> rerun;
        std::cout << "Please hit Enter when ready.\n";
        system("pause");
    }
    while( rerun == 'Y');

    return 0;
}