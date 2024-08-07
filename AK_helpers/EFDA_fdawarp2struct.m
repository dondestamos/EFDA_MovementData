function warpstruct = fdawarp2struct(warpobj)
warpstruct = struct;
fldns = fieldnames(warpobj);
for ifldn = 1:length(fldns)
    fldn = fldns{ifldn};
    warpstruct.(fldn) = warpobj.(fldn);
end
end