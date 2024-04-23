%Ivan Arias
%2017/10/20
%Convert data to NetCDF
%Run in gannett

%To convert RAW Vaisala Data into netCDF

DataPath_month = '/net/denali/storage2/radar2/tmp/RELAMPAGO/RAW/2019/01/';
OutputDirectory = '/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/';
RadxConvertPath = '/usr/local/bin/RadxConvert';
Directory_days = dir(DataPath_month);

for day = 29:length(Directory_days)
    DataPath = [DataPath_month '/' Directory_days(day).name];
    cd(DataPath)
    Directory = dir('COL*');
    for i = 1:length(Directory)
        str = [num2str(i) '/' num2str(length(Directory))];
        disp(str)
        try
            filename = Directory(i).name;
            cmd = [RadxConvertPath ' -cfradial ' ' -f ' filename ' -outdir ' OutputDirectory];
            status = system(cmd);
            if status == 0
               str = ['file ' Directory(i).name ' done'];
               disp(str);
            else
               str = ['Unable to process ' Directory(i).name];
               disp(str)
            end
        catch 
            str = ['An error occurred while processign ' Directory(i).name];
            disp(str);
        end
    end
end
