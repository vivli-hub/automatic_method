function run_filt_WHO(wreg, flag_transf, flag_plot, sp, ftype, fdomain, window)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%window         'hann', 'halfhann', or 'rect'
%%%%%%%%%%%%%%%%%%%%%%f%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%load list countries
% load('file_list_countries_short.mat','list_countries');
list_countries   = {
'Afghanistan'
'Albania'
'Algeria'
'American Samoa'
'Andorra'
'Angola'
'Anguilla'
'Antigua and Barbuda'
'Argentina'
'Armenia'
'Aruba'
'Australia'
'Austria'
'Azerbaijan'
'Bahamas'
'Bahrain'
'Bangladesh'
'Barbados'
'Belarus'
'Belgium'
'Belize'
'Benin'
'Bermuda'
'Bhutan'
'Bolivia (Plurinational State of)'
'Bonaire'
'Bosnia and Herzegovina'
'Botswana'
'Brazil'
'British Virgin Islands'
'Brunei Darussalam'
'Bulgaria'
'Burkina Faso'
'Burundi'
'Cabo Verde'
'Cambodia'
'Cameroon'
'Canada'
'Cayman Islands'
'Central African Republic'
'Chad'
'Chile'
'China'
'Colombia'
'Comoros'
'Congo'
'Cook Islands'
'Costa Rica'
'Ctedvoire'
'Croatia'
'Cuba'
'Curaao'
'Cyprus'
'Czechia'
'Democratic Peoples Republic of Korea'
'Democratic Republic of the Congo'
'Denmark'
'Djibouti'
'Dominica'
'Dominican Republic'
'Ecuador'
'Egypt'
'El Salvador'
'Equatorial Guinea'
'Eritrea'
'Estonia'
'Eswatini'
'Ethiopia'
'Falkland Islands (Malvinas)'
'Faroe Islands'
'Fiji'
'Finland'
'France'
'French Guiana'
'French Polynesia'
'Gabon'
'Gambia'
'Georgia'
'Germany'
'Ghana'
'Gibraltar'
'Greece'
'Greenland'
'Grenada'
'Guadeloupe'
'Guam'
'Guatemala'
'Guernsey'
'Guinea'
'Guinea-Bissau'
'Guyana'
'Haiti'
'Holy See'
'Honduras'
'Hungary'
'Iceland'
'India'
'Indonesia'
'Iran (Islamic Republic of)'
'Iraq'
'Ireland'
'Isle of Man'
'Israel'
'Italy'
'Jamaica'
'Japan'
'Jersey'
'Jordan'
'Kazakhstan'
'Kenya'
'Kiribati'
'Kosovo'
'Kuwait'
'Kyrgyzstan'
'Lao Peoples Democratic Republic'
'Latvia'
'Lebanon'
'Lesotho'
'Liberia'
'Libya'
'Liechtenstein'
'Lithuania'
'Luxembourg'
'Madagascar'
'Malawi'
'Malaysia'
'Maldives'
'Mali'
'Malta'
'Marshall Islands'
'Martinique'
'Mauritania'
'Mauritius'
'Mayotte'
'Mexico'
'Micronesia (Federated States of)'
'Monaco'
'Mongolia'
'Montenegro'
'Montserrat'
'Morocco'
'Mozambique'
'Myanmar'
'Namibia'
'Nauru'
'Nepal'
'Netherlands'
'New Caledonia'
'New Zealand'
'Nicaragua'
'Niger'
'Nigeria'
'Niue'
'North Macedonia'
'Northern Mariana Islands (Commonwealth of the)'
'Norway'
'occupied Palestinian territory including east Jerusalem'
'Oman'
'Other'
'Pakistan'
'Palau'
'Panama'
'Papua New Guinea'
'Paraguay'
'Peru'
'Philippines'
'Pitcairn Islands'
'Poland'
'Portugal'
'Puerto Rico'
'Qatar'
'Republic of Korea'
'Republic of Moldova'
'Runion'
'Romania'
'Russian Federation'
'Rwanda'
'Saba'
'Saint Barthlemy'
'Saint Helena Ascension and Tristan da Cunha'
'Saint Kitts and Nevis'
'Saint Lucia'
'Saint Martin'
'Saint Pierre and Miquelon'
'Saint Vincent and the Grenadines'
'Samoa'
'San Marino'
'Sao Tome and Principe'
'Saudi Arabia'
'Senegal'
'Serbia'
'Seychelles'
'Sierra Leone'
'Singapore'
'Sint Eustatius'
'Sint Maarten'
'Slovakia'
'Slovenia'
'Solomon Islands'
'Somalia'
'South Africa'
'South Sudan'
'Spain'
'Sri Lanka'
'Sudan'
'Suriname'
'Sweden'
'Switzerland'
'Syrian Arab Republic'
'Tajikistan'
'Thailand'
'The United Kingdom'
'Timor-Leste'
'Togo'
'Tokelau'
'Tonga'
'Trinidad and Tobago'
'Tunisia'
'Trkiye'
'Turkmenistan'
'Turks and Caicos Islands'
'Tuvalu'
'Uganda'
'Ukraine'
'United Arab Emirates'
'United Republic of Tanzania'
'United States of America'
'United States Virgin Islands'
'Uruguay'
'Uzbekistan'
'Vanuatu'
'Venezuela (Bolivarian Republic of)'
'Viet Nam'
'Wallis and Futuna'
'Yemen'
'Zambia'
'Zimbabwe'
};
%folder data
folder_covid    = './WHO_Data/';

if strcmp(sp, 'week')
    T0 = 52;
    n = 4;
    delta = 1/10;
elseif strcmp(sp, 'day')
    T0 = 7;
    n = 6;
    delta = 1/7;
else
end

%Initialization
cc = delta*T0;
Wn = [1/(T0+cc) 1/(T0-cc)]; %for 'stop'
fn = 1/(T0/2 + cc);         %for 'high'
ah = [];
seven = T0./(1:5);
BICo = zeros(1,length(seven));
BICf = BICo;
SCo = BICo;
SCf = BICf;
noc = 0; %no. of countries

%Files for saving the results
fina = strcat('Results_', wreg, '_tr', num2str(flag_transf), '_', sp, '_', ftype, '_', fdomain, '_', window);
t_before={}; 
t_after={};
save(strcat(fina,'.mat'), 't_before', 't_after', 'seven');
fid = fopen(strcat(fina,'.txt'),'w');


for i=1:length(list_countries)

fname = strcat(folder_covid,list_countries{i},'.mat');
load(fname,'y','whoreg');

if strcmp(wreg,whoreg) %region   
    noc = noc+1;   
%Vector y can have negative entries and NaN!!!
ytr         = y;
ytr(ytr<0)  = 0;
ytr(isnan(ytr))  = 0;
if sum(abs(ytr)) == 0 %all entries are equal to zero
    fprintf(fid, '%s \t %s \t %s\n', list_countries{i}, whoreg, 'No data');
else
    
if strcmp(sp,'week')
    [ytr, err] = data_week(ytr); 
       if err==1 
          return 
       end
end

%Apply transform
ytr = transf_data(ytr, flag_transf);

%Initialization
[A,vec_om,vec_f] = init(ytr);

%Spectral estimation (before filtering)
[Phi_WP, freq_BIC, freq_NML] = get_spect(ytr,vec_f,vec_om,A,window);
if numel(union(freq_BIC,freq_NML))==0
    fprintf(fid, '%s \t %s \t %s\n', list_countries{i}, whoreg, 'No periodicities');
else
    [BICo, SCo] = save_res('before_filt',seven, fina, fid, list_countries{i}, whoreg, freq_BIC, freq_NML, T0, delta, BICo, SCo);
end

if flag_plot==1
    label = strcat(list_countries{i}, ' (before filtering)');
    fig = plot_periodogram(1, 1, Phi_WP(1:length(vec_f)), vec_f, freq_BIC, freq_NML, label, sp);
end

if sum( (freq_BIC>=min(Wn)).*(freq_BIC<=max(Wn)) ) + sum( (freq_NML>=min(Wn)).*(freq_NML<=max(Wn) ) )  >0  

    %Design filter
    if length(ah) ~= length(ytr) 
        [ah,b,a] = getfilt(length(ytr),Wn,fn,n,ftype);
    end
   
    %Apply filter
    %Freq. domain
     te = ah.*fft(ytr); 
     if max( abs( te - conj( te([1,end:-1:2]) ) ) )>0
         fprintf('Error!\n');
         return
     end
     yf = ifft(te,'symmetric');
     %Time domain
     yf2 = filter(b,a,ytr(end:-1:1));
     yf2 = filter(b,a,yf2(end:-1:1));
    
     if flag_plot==1
        figure(2)
        plot(1:length(ytr),yf2,'r-',1:length(ytr),yf,'b-',1:length(ytr),ytr,'k:');
        legend('Time domain filter','Freq. domain filter','Original time series');
     end
     
    %Spectral estimation (after filtering)
    if strcmp(fdomain,'freq')
        [Phi_WP_f, freq_BIC_f, freq_NML_f] = get_spect(yf,vec_f,vec_om,A,window);
    elseif strcmp(fdomain,'time')
        [Phi_WP_f, freq_BIC_f, freq_NML_f] = get_spect(yf2,vec_f,vec_om,A,window);
    else
        fprintf('Error: Incorrect value fdomain');
        return;
    end
    [BICf, SCf] = save_res('after_filt', seven, fina, fid, list_countries{i}, whoreg, freq_BIC_f, freq_NML_f, T0, delta, BICf, SCf);
    if flag_plot==1
        label = strcat(list_countries{i}, ' (after filtering)');
        plot_periodogram(1, 2, Phi_WP_f(1:length(vec_f)), vec_f, freq_BIC_f, freq_NML_f, label, sp);
    end

end %if max

if flag_plot==1
    pause
    close(fig);
end

end % if sum(abs(ytr))

end %if strcmp(wreg,whoreg)
    
end %for list_countries

save(strcat(fina,'.mat'),'noc','BICo','BICf','SCo','SCf', '-append');

fclose(fid);

close all

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ytr = transf_data(y, flag_transf)
    
    if flag_transf==0
    elseif flag_transf==1
        y = y(2:end)-y(1:end-1);
        y(y==0) = 1;
        y = sign(y).*log(abs(y));
        y = [0; y];
    elseif flag_transf==2
        y(y==0) = 1;
        y = 100*( log( y(2:end) ) - log( y(1:end-1) ) );
        y = [0; y];
     elseif flag_transf==3
         y = y(2:end)-y(1:end-1);
         y = [0; y];
    else
        fprintf('Error: Incorrect value for flag_transf.\n');
        return;
    end
    
    ytr = y;
    
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,vec_om,vec_f] = init(ytr)

            N = length(ytr);
            T = (1:N);
            T = T(:);
            vec_om = 2*pi*(0:1/N:1/2);
            vec_f = vec_om/(2*pi);
            K = length(vec_om);
            C = cos(kron(vec_om,T));
            S = sin(kron(vec_om,T));
            A = zeros(N,2*K);   
            A(:,1:2:end) = C;
            A(:,2:2:end) = S;
            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi_WP, freq_BIC, freq_NML] = get_spect(ytr,vec_f,vec_om,A,window)

            L = length(ytr);

            if strcmp(window,'hann')
                win = hanning(L); 
            elseif strcmp(window,'halfhann')
                win = halfhann_w(L);
            elseif strcmp(window,'rect')
                win = ones(1,L);
            else
                fprintf('Error in get_spect function\n');
            end

            ytr = detrend(ytr);
            
            [Phi_WP,f] = periodogram(ytr,win,(0:1/L:1/2)',1);
            
            if f~=vec_f
                fprintf('Error in get_spect function\n');
            end
            Phi_WP(Phi_WP<=0) =  min( Phi_WP(Phi_WP>0) )/10;

            Mmax = 20;
            [~,ind] = sort(Phi_WP(1:length(vec_f)),'descend');
            
            rho = 4;
            M_BIC = comp_BIC(Mmax,ind,rho,ytr,A);
            freq_BIC = vec_f(ind(1:M_BIC));
            
            [M_NML,~] = comp_NML(Mmax,ind,ytr,A,vec_om);
            freq_NML = vec_f(ind(1:M_NML));
            
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = halfhann_w(M)
    temp = (0:M-1);
    w = 0.5 + 0.5*cos((pi.*temp)./M);
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M_BIC = comp_BIC(Mmax,ind,rho,y,A)

N = length(y);

BIC = zeros(Mmax+1,1);
err = y;
RSS =  sum(err.^2);
BIC(1) = N*log(RSS);
Ak = [];
for M=1:Mmax
    k = ind(M);
    Ak = [Ak A(:,2*k-1:2*k)];
    if cond(Ak)<10^3
        Teta = Ak\y;
        err = y - Ak*Teta;
        RSS = sum(err.^2);
        BIC(M+1) = N*log(RSS) + rho*M*log(N);
    else
        BIC(M+1) = Inf;
    end
end

[~,w]=min(BIC);
M_BIC = w-1;

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M_NML,diff_pen] = comp_NML(Mmax,ind,y,A,vec_om)

PEN = zeros(1,Mmax+1);

N = length(y); %no. of samples

det_Q2 = N^5;
G1 = [1/3 1/2; 1/2 1];
det_G1 = det(G1);

NML = zeros(Mmax+1,1);

err = y;
RSS =  sum(err.^2);
%compute det(FIM)
s2 = RSS/N;
det_FIM_s2 = N/2/s2^2;
det_FIM = det_FIM_s2;
%compute NML
NML(1) = (N/2)*log(RSS) + log(det_FIM)/2 + log(s2+N^(-1/4));
PEN(1) = 2*( NML(1) - (N/2)*log(RSS) );

om_vec = [];
Ak = [];
for M=1:Mmax
 
    k = ind(M);
    Ak = [Ak A(:,2*k-1:2*k)];
    if cond(Ak)<10^3
        Teta = Ak\y;
        err = y - Ak*Teta;
        RSS = sum(err.^2);
    
        [alpha,phi]=change_param(Teta);
   
        alpha_vec = alpha;
        phi_vec = phi;
        om_vec = [om_vec; vec_om(k)];
    
        %FIM s2
        s2 = RSS/N;
        det_FIM_s2 = N/2/s2^2;
        %FIM alpha,omega,phi
        snr_vec = (alpha_vec.^2)/2/s2;
        det_G = (prod(snr_vec)^3)*(prod(alpha_vec)^(-2))*(det_G1^M);
        %FIM
        det_FIM = det_FIM_s2*(det_Q2^M)*det_G;
    
        temp_vec = [alpha_vec; phi_vec; om_vec; s2];
        temp = sum( log( abs(temp_vec)+N^(-1/4) ) );
        NML(M+1) = (N/2)*log(RSS) + log(det_FIM)/2 + temp;
        PEN(M+1) = 2*( NML(M+1) - (N/2)*log(RSS) );
    else
        NML(M+1) = Inf;
        PEN(M+1) = Inf;
    end
end

[~,w]=min(NML);
if w<length(NML)
    diff_pen = PEN(w+1)-PEN(w);
else
    diff_pen = 0;
end
M_NML = w-1;

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alphaT,phiT]=change_param(TetaT)

alphaT = [];
phiT = [];

for i=1:2:size(TetaT,1)/2
    
    Teta = TetaT(2*i-1:2*i,:);

    ak = Teta(1,:);
    bk = Teta(2,:);

    phi = -atan(bk./ak);
    ff = find(sin(phi)==0);
    rff = setdiff(1:length(phi),ff);
    alpha = zeros(size(phi));
    alpha(rff) = -bk(rff)./sin(phi(rff));
    alpha(ff) = ak(ff)./cos(phi(ff));
    % alpha is positive
    ff = find(alpha<0);
    alpha(ff) = -alpha(ff);
    ss=(sign(phi(ff))>=0)-1/2;
    phi(ff) = phi(ff)-2*ss*pi;

    alphaT = [alphaT; alpha(:)];
    phiT = [phiT; phi(:)];
    
end 

end %function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig = plot_periodogram(fig, panel, Periodogram, Freq, freq_BIC, freq_SC, label, sp)

fig = figure(fig);
subplot(1,2,panel);

Periodogram(Periodogram<10^(-15)) = min( Periodogram(Periodogram>10^(-15)) )/10;
Periodogram1 = Periodogram;
Freq1 = Freq;

%periodogram
plot(Freq1,10*log10(Periodogram1),'k-');
xlim([0 0.1]);
hold on

%BIC estimates
if numel(freq_BIC)>0
    pos = get_pos(Freq1, freq_BIC);
    val = 10*log10(Periodogram1(pos));
    plot(Freq1(pos),val,'ro');
    legend('BIC');
    hold on
end

%SC estimates
if numel(freq_SC)>0
    pos = get_pos(Freq1, freq_SC);
    val = 10*log10(Periodogram1(pos));
    plot(Freq1(pos),val,'b*');
    if numel(freq_BIC)>0
        legend('Periodogram','BIC','SC','Location','southeast');
    else
        legend('Periodogram','SC','Location','southeast');
    end
    hold on
end

xlim([0 1/2]);
xl = strcat('Frequency [cycles/', sp, ']');
xlabel(xl);
% xlabel('Frequency [cycles/week]');
ylabel('Amplitude [dB]');
title(label);

hold off

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = get_pos(Freq1, freq)

lf = length(freq);
pos = zeros(lf,1);
for i=1:lf
    diff = Freq1-freq(i);
    [~,p] = min(abs(diff));
    pos(i) = p;
end


end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ah,b,a] = getfilt(N,Wn,fn,n,ftype) 

Rp = 0.01;
Rs = 40;

if strcmp(ftype,'stop')
    % Transfer function design
    [b,a] = ellip(n,Rp,Rs,Wn./(1/2),ftype);
elseif strcmp(ftype,'high')
    % Transfer function design
    [b,a] = ellip(n,Rp,Rs,fn/(1/2),ftype);
else
    fprintf('Error: Incorrect filter type\n');
    return;
end

ah = fft(b,N)./fft(a,N);
ah = ah(:);

end %function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ytr, err] = data_week(ytr)

err = 0;
rest = rem(length(ytr),7);
ymat = reshape( ytr(rest+1:end), [7,length(ytr(rest+1:end))/7] );
yrow = sum(ymat,1);
ytr = yrow(:);
if length(ytr)>512
    fprintf('Error\n');
    ytr = []; err = 1;
    return
else
    len = length(ytr)-256;
    ytr1 = ytr(1:len);
    ytr2 = ytr(len+1:end);
    if sum(abs(ytr1)~=0)
        fprintf('Error\n');
        ytr = []; err = 1;
    else
        ytr = ytr2;
    end
end %if length
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BIC, SC] = save_res(label_filt, seven, fina, fid, country, whoreg, freq_BIC, freq_NML, T0, delta, BIC, SC)

t{1} = {country};
for i=1:length(seven)
    t{i+1} = {};
end

fprintf(fid, '%s \t %s', country, whoreg);

Cm = unique(union(1./freq_BIC,1./freq_NML));
Cm = sort(Cm,'descend');

temp = [];
for i=1:length(Cm)
    [mm,w] = min( abs(seven - Cm(i)) );
     if mm<( delta*T0 )
         t{w+1} = [t{w+1} Cm(i)];
         if ~ismember(seven(w),temp)
            temp = [temp seven(w)];
            fprintf(fid, '\t %i*%s', T0, strtrim(rats(seven(w)/T0)));
            
            fprintf(fid,'[');
            if ismember(Cm(i),1./freq_BIC) 
                fprintf(fid,'%s','B');
                BIC(w) = BIC(w)+1;
            end
            if ismember(Cm(i),1./freq_NML)
                fprintf(fid,'%s','S');
                SC(w) = SC(w)+1;
            end
            fprintf(fid,']');
         end
     end
end

fprintf(fid, '\n');

if strcmp(label_filt,'before_filt')
    load(strcat(fina,'.mat'),'t_before');
    t_before = [t_before; t];
    save(strcat(fina,'.mat'),'t_before','-append');
elseif strcmp(label_filt,'after_filt')
    load(strcat(fina,'.mat'),'t_after');
    t_after = [t_after; t];
    save(strcat(fina,'.mat'),'t_after','-append');
else
end
         
end %function





