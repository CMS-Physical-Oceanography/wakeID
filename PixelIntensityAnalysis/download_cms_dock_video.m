  % Define the URL for the desired day: yyyy/mm/dd/                         
    url = 'https://stage-ams.srv.axds.co/archive/mp4/uncw/cms_dock_south/2022/10/15/';               
                                                                              
    % read the list of files within this folder                               
    urlData = webread([url]);                                                 
    % get rid of all the html code                                            
    fileList = regexp(urlData,'(?<=<pre>).*(?=</pre>)','match');              
    % find the start and end of each filename (start=cms, end=mp4)            
    cmsInds = regexp(fileList{1},'cms');                                      
    mp4Inds = regexp(fileList{1},'mp4');                                      
    % loop over the start/end indices to get filenames                        
    for ii = 1:length(cmsInds); files{ii} =                                   
  fileList{1}(cmsInds(ii):mp4Inds(ii)); end                                   
                                                                              
    % download the first video to your Downloads folder                       
    fileOut = ['~/Downloads/',files{1}];                                      
    websave(fileOut, [url, files{1}]);  
