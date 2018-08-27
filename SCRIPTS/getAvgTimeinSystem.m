function avgTimeinSystem = getAvgTimeinSystem(userList,simSlotsPerSec)
    TinSlots = 0;
    nOfUsers = length(userList);
    for i = 1 : nOfUsers
        if(userList(i).state == 2)
            TinSlots = TinSlots +  userList(i).exitTime - userList(i).entryTime;
        end
    end
    avgTimeinSystem = TinSlots / (nOfUsers * simSlotsPerSec);
end