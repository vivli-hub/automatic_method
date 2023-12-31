function read_countries_latitude

folder_covid    = './WHO_Data/';
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

%Read from csv to mat
if ~isfile('T_WHO_FILE.mat')
    csvfile = strcat(folder_covid,'WHO_data.csv');
    fprintf('The data are from the file %s\n', csvfile);
    TT = readtable(csvfile);
    T = table2cell(TT); 
    save('T_WHO_FILE.mat','T');
else
    fprintf('The data are from the file T_WHO_file.mat\n');
    load('T_WHO_file.mat','T');
end

%Read
maxd = 0;
for ind=1:length(list_countries)

        c   = list_countries(ind);
        ff  = find(strcmp(T(:,4),c)); %name of the country is on the 3rd col.
        if isempty(ff)
            fprintf('Data for %c does not exist\n',char(c));
        else
            cd = T(ff,6); %no. of daily cases is on the 5th col.
            d = Inf(length(ff),1);
            for i=1:length(ff)
                d(i) = cd{i,1};
            end
            
            d(d<0) = 0;         %Replace negative numbers with zeros
            d(isnan(d)) = 0;    %Replace NaN with zeros
            
            if length(d)>maxd
                maxd=2^ceil( log2(length(d)) );
            end
            
            %Write (one mat-file for each country)
            y = d;
            if ~isempty(y)
                if length(y)<maxd
                    y = [zeros(maxd-length(y),1); y];
                end
            end
            
            whoreg = char( T(ff(1),10) ); %read WHO reg. from col. 4
            fname = strcat(folder_covid,list_countries{ind},'.mat');
            save(fname,'y','whoreg');
        end
end %ind

end %function

