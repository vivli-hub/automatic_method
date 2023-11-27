
# automatic_method
An automatic method for the identification of cycles in Covid-19 time series data

## Background
In May 2023, the World Health Organization (WHO) declared the end of the global health emergency related to COVID-19. Since January 2020, when WHO classified Covid-19 as a public health emergency of international concern, a substantial amount of epidemiological data has been gathered. Throughout this timeframe, numerous attempts have been made to analyze the data. However, many of these investigations were constrained due to a lack of sufficient measurements available at the time of analysis. There are various aspects that warrant exploration, with one of them focusing on identifying "waves" in the pandemic. This involves uncovering cycles in the collected data for each country/territory/area (CTA). Given the absence of a comprehensive study for all 236 CTAs with available time series data, we undertake this analysis in the current study.

### Abstract
All previous methods for the identification of cycles in Covid-19 daily and weekly data involve a subjective interpretation of the results. This poses difficulties for researchers interested in conducting a comprehensive study which analyzes the presence of the cycles for each CTA. We have designed an algorithm that detects automatically the fundamental period T0 and the harmonics T0/2,...,T0/5, where T0=7 days for daily data and T0=52 weeks for weekly data. We have tested the new algorithm by applying it to the time series from 236 CTA's, where WHO collected the Covid-19 data. The detection results we have obtained confirm the findings previously reported by other researchers.
In this talk, we will present all the details of our algorithm and comment on the results obtained in experiments with Covid-19 time series data. We will also discuss a proposal for evaluating the dissimilarity between the time series collected for two different CTA’s.

## Data
In our experiments, we use the daily counts of COVID-19 cases that are publicly available on the website of WHO: [Data](https://covid19.who.int/data). More precisely, we consider $x_1,\ldots,x_n$ to be the time series that represents the number of daily cases from 3 January 2020 to 3 May 2023, for a particular CTA. 

## automatic_method
- [read_countries_WHO](#read_countries_WHO)
- [run_filt_WHO](#run_filt_WHO)
- [run_filt_fullsearch](#run_filt_fullsearch)
- [read_countries_latitude](#read_countries_latitude)
- [run_filt_latitude](#run_filt_latitude)

### read_countries_WHO

Data can be downloaded from: [Data](https://covid19.who.int/WHO-COVID-19-global-data.csv). In the folder WHO/files/Data there are one data file WHO_data.csv. This file contains the data that I have used in my experiments. It is good to use the same data when you will try to understand the code that I have wrote. This data file is in the WHO/files/Data. Function read_countries_WHO.m reads the data from WHO_data.csv and creates a Mat file for each country in the folder WHO/files/Matlab_files/WHO_Data. I have created these files for all countries (There are 236 countries in my experiments.)

### run_filt_WHO

Once you have these Mat files, you can run the main function run_filt_WHO.m. You need to give as input a flag which indicates how the data are transformed (have a look at the comments in the function to understand how to select this flag)
For example, you can call the function like this:
run_filt(wreg, flag_transf, flag_plot, sp, ftype, fdomain)
%Input:
%wreg:          WHO region for which the data is analysed    
%flag_transf:   flag for data transformation
%               0=no transf. is applied to the data
%               1=transf. from GP Nason, Scientific Reports, 2020
%               2=transf. which is usually applied to daily stock markets indices at closing time
%               3=first order differences (computed in time domain)
%flag_plot:     0=no plots, 1=plots for each country
%sp             'day'=daily data, 'week'=weekly data
%ftype          'stop'=stop band filter, 'high'=high pass filter
%fdomain        'time'=filtering in time domain, 'freq'=filtering in freq. domain
run_filt_WHO(wreg, flag_transf, flag_plot, sp, ftype, fdomain)
Note that the WHO region is now an input parameter (wreg).

run_filt_WHO produces two different files. Both of them have the same name, which is given by fina = strcat('Results_', wreg, '_tr', num2str(flag_transf), '_', sp,
'_', ftype, '_', fdomain); However, one file has the extension .txt and the other file has the extension .mat. The .txt file is very similar to the file produced by the previoud version of the code. The .mat file contains the following variables:

seven: T0./(1:5);
t_before: is a matrix of size n*5, where n represents the number of countries in a particular region that can be detected with at least one periodogram or harmonics before applying the filter. The 5 columns correspond to the detection results from T0 to T0/5, where the numbers represent the periods associated with frequencies detected by BIC or SC. It's worth noting that the n in this context does not represent the number of countries contained within a particular region, as there are some countries that cannot be detected with any periodogram or harmonics；
t_after: Most of the content in this also has the same meaning as "t_before." The difference between them is that "t_after" stores the data conclusions after applying the filter. Only countries that were detected with a fundamental periodogram before applying the filter will apply the filter. 

### run_filt_fullsearch

Once you have these Mat files, you can run the main function run_filt_fullsearch.m. You need to give as input a flag which indicates how the data are transformed (have a look at the comments in the function to understand how to select this flag)

%Input:
%wreg:          WHO region for which the data is analysed    
%flag_transf:   flag for data transformation
%               0=no transf. is applied to the data
%               1=transf. from GP Nason, Scientific Reports, 2020
%               2=transf. which is usually applied to daily stock markets indices at closing time
%               3=first order differences (computed in time domain)
%
%sp             'day'=daily data, 'week'=weekly data
%ftype          'stop'=stop band filter, 'high'=high pass filter
%fdomain        'time'=filtering in time domain, 'freq'=filtering in freq. domain


The output generated by this function is very similar to that of run_filt_WHO.m. Both are algorithms for detecting the fundamental periodogram and harmonics, with slight differences in their processes. The specific distinctions are explained in the PDF document.

### read_countries_latitude

Function read_countries_latitude.m reads the data from WHO_data.csv and creates a Mat file for each country in the folder Latitude/files/Matlab_files/WHO_Data
It's worth noting that the WHO website does not provide specific information about the latitude position for each country. Therefore, we need to carry out some manual labeling. I have created these files for all countries (There are 236 countries in my experiments.) we extend our analysis by categorizing the CTA’s according to the latitudinal position of the capital city. This relies on the fact that the capital cities are typically in central locations within countries and have the highest economics activities and population densities.


### run_filt_latitude

Once you have these Mat files, you can run the main function run_filt_latitude.m. You need to give as input a flag which indicates how the data are transformed (have a look at the comments in the function to understand how to select this flag)
For example, you can call the function like this:
run_filt(latitude, flag_transf, flag_plot, sp, ftype, fdomain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%latitude:          Latitude group for which the data is analysed    
%flag_transf:   flag for data transformation
%               0=no transf. is applied to the data
%               1=transf. from GP Nason, Scientific Reports, 2020
%               2=transf. which is usually applied to daily stock markets indices at closing time
%               3=first order differences (computed in time domain)
%flag_plot:     0=no plots, 1=plots for each country
%sp             'day'=daily data, 'week'=weekly data
%ftype          'stop'=stop band filter, 'high'=high pass filter
%fdomain        'time'=filtering in time domain, 'freq'=filtering in freq. domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_filt_latitude produces two different files. Both of them have the same name, which is given by fina = strcat('Results_', wreg, '_tr', num2str(flag_transf), '_', sp,
'_', ftype, '_', fdomain);
However, one file has the extension .txt and the other file has the extension .mat. The .txt file is very similar to the file produced by the previous version of the code.
The .mat file contains the following variables:

noc = total number of countries in the particular latitude group
seven = seven = T0./(1:5);
BICo has five entries:
first entry: total no. of the countries in the latitude group for which BIC has found the periodicity T0 in the data, before filtering
second entry: total no. of the countries in the latitude group for which BIC has found the periodicity T0/2 in the data, before filtering
fifth entry: total no. of the countries in the latitude group for which BIC has found the periodicity T0/5 in the data, before filtering
BICf   = similar to BICo, except that the statistics are for the filtered data
SCo = similar to BICo, except that SC is used instead of BIC
SCf = similar to BICf, except that SC is used instead of BIC


The full report can be found at [Link]().
