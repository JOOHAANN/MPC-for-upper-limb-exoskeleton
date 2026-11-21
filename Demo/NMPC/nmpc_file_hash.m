function value=nmpc_file_hash(path)
%NMPC_FILE_HASH Identify the exact model/code revision used by cached runs.
file=fopen(path,'rb');assert(file>=0);cleanup=onCleanup(@()fclose(file));
engine=java.security.MessageDigest.getInstance('SHA-256');
engine.update(fread(file,inf,'*uint8'));
value=lower(reshape(dec2hex(typecast(engine.digest(),'uint8'),2)',1,[]));
end
