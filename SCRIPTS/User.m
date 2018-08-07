classdef User
    % User stores details of a DASH user 
    %  We store all simulation results for a user in an object of class User 
    
    properties
        id;
        class;
        videoLength; % in second (should be exponentially distributed)
        totalNbrOfSegments;
        bufferLevel=0;
        currentSegment=0;
        slotsUsed = 0;% an attribute that stores the number of simulation slots used, 
        %ideal use case is to find slots used for each segment then update it to zero for the next one.
        prevQuality=0;
        currentQuality;
        qualitySwitches=0;
        avgQuality = 0;
        propOfSegmentDownloaded=0;
        state=0; % states are "requesting" == 0, "streaming" == 1, "outOfSystem" == 2, "balk" == 3
        isPlaying=0;
        starv=0;
        starvedSegment = -1;
        visitsBadState = 0;
        trace=trace_struct;
        usersInSystem; % Number of existing users seen by arriving user, can be associated with its prefetching delay too
        stateAtPrefetchStart; % this is the state of the system when this user started prefetching 
        inPrefetchStart = 1; % if user is in prefetch beginnning then this is 1 else 0
        updateDelayMatrixFlag = 0; 
        prefetchTime = -1;
        entryTime = -1;
        exitTime = -1;
        prefFlag = 1;
    end    
    methods     
        
    end    
end
