for i = 1 : 270
    if abs(D.eventtime{i}(3) - D.eventtime{i}(2)) > 1
        repeat2(i,1) = 1;
    else
        repeat2(i,1) = 0;
    end
end