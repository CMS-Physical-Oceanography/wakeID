clear all
close all

numdays = 10;
daystart = 1;

% read the list of files within this folder                               
% urlData = webread([url]);                                                 
% get rid of all the html code                                            
% fileList = regexp(urlData,'(?<=<pre>).*(?=</pre>)','match');              
% find the start and end of each filename (start=cms, end=mp4)            
% cmsInds = regexp(fileList{1},'cms');                                      
% mp4Inds = regexp(fileList{1},'mp4')+2;                                      
% loop over the start/end indices to get filenames                        
% for ii = 1:length(cmsInds);
   % files{ii} = fileList{1}(cmsInds(ii):mp4Inds(ii));
    
% $$$     % download the first video to your Downloads folder                       
% $$$     fileOut = ['~/Downloads/',files{ii}];                                      
% $$$     websave(fileOut, [url, files{ii}]);  
% end                                   

baseURL = 'https://stage-ams.srv.axds.co/archive/mp4/uncw/cms_dock_south/2022/10/';

for k = 1:numdays

    dayofmonth = daystart + k - 1;
    dayStr = sprintf('%02d', dayofmonth);

    downloadsFolder = fullfile('E:\Boatvideos','October',dayStr);

    if ~exist(downloadsFolder,'dir')
        mkdir(downloadsFolder)
    end

    url = [baseURL, dayStr, '/'];

    urlData = webread(url);

    fileList = regexp(urlData, '(?<=<pre>).*(?=</pre>)', 'match');

    files = regexp(fileList{1}, 'cms.*?\.mp4', 'match');

    for ii = 1:length(files)

        fileOut = fullfile(downloadsFolder, files{ii});
        fileURL = [url, files{ii}];

        fprintf('Downloading %s\n', fileURL)

        websave(fileOut, fileURL);

    end
end