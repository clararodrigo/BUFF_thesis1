% Bingxue for data serialization and saving
function saveIQData_simple(data, filename, pixelMap, UserSet)

    S = whos('data');
    splitPart = ceil(S.bytes/2^32);

%     n = ndims(data);
%     for i = 2 : splitPart
%         part = ceil((i-1)*size(data,n)/splitPart+1):ceil(i*size(data,n)/splitPart);
%         eval(['data',num2str(i),'=data(:,:,part);']);
%     end
%     part = 1:ceil(1*size(data,n)/splitPart);
%     data = data(:,:,:,part);

%     data = serialize(data);
    data = hlp_serialize(data);
%     if splitPart ~= 1
%         for i = 2 : splitPart
%             eval(['data',num2str(i),'=hlp_serialize(data',num2str(i),');']);
%         end

%         string = ['savefast([savedir,name,''_BF_'' ImagParam.imageType ''_'' ImagParam.returnType ''.mat''],''data'',''UserSet'',''ImagParam'',''pixelMap'',''splitPart'''];
%         for i = 2:splitPart
%             var=  ['''data',num2str(i),''''];
%             string = [string,',',var];
%         end
%         string=[string,')'];
%         eval(string);
%     else
        savefast(filename,'data','UserSet','pixelMap');
        disp('done saving');
%     end

    clear splitPart;

end