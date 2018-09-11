function [user, update] = DASH_KR(paramsDASH, segPerSlot, user, time, C_time, bThres, idThres, ScalarSystemState)
% The DASH code 
% This function simulates the DASH code 
% It is run every simulation slot and in each sim. slot it
    % downloads some video  
    % plays out some video (if user is playing)
    % prefetches some video (if user is prefetching)
% In each slot, we download video till there is some capacity left.
% Also the player plays out video or prefetches video (depending on playout buffer status)
% If player sees starvation, we record a stravation event and go to prefetching mode.
% We make a decision on the video bitrate for each segment based on the playout buffer value.
% If user finishes service, it stores all its results in its corresponding object.
downloaded = 0;
%totalC = C_time; % total bits available for this user in current sim slot
epsilon = 1e-5;
propSlot = 0;
update = 0;
% if user.currentSegment <= user.totalNbrOfSegments
if (user.state == 0)
    if (user.currentSegment == user.totalNbrOfSegments)
        user.state = 2;%out of system
        user.avgQuality = user.avgQuality / user.totalNbrOfSegments;
        user.exitTime = time;% stores the slot in which the user exits
        if (user.prefetchTime <= 0) % Video was delivered completely in prefetching.
            user.prefetchTime = time; % check this
        end
        update = -1;
        return;
    end
    user.currentSegment = user.currentSegment + 1;
    user.propOfSegmentDownloaded = 0;%set to zero for this current segment
    
    % Make decision on the next quality
    % (buffer based decision)  STARTS HERE
    %printf('%duserBuffer \n',user.bufferLevel, paramsDASH.Bmin);
    %printf('%dBmin \n', user.bufferLevel, paramsDASH.Bmin);
    if (user.bufferLevel <= paramsDASH.Bmin)% paramsDASH.Bmin is 4
        user.currentQuality = paramsDASH.l(1);%the minimum bit rate
    else
        if (length(bThres) > 1) % DISCRETE quality levels
            idx = max(idThres(bThres <= user.bufferLevel));
            user.currentQuality = paramsDASH.l(idx);
        else                    % CTX quality levels
            user.currentQuality = paramsDASH.rate * user.bufferLevel + paramsDASH.c;
            user.currentQuality = min(paramsDASH.l(end), user.currentQuality);
        end
    end
    if(user.prefFlag) % throughout pref give lowest quality.
        user.currentQuality =  paramsDASH.l(1);
    end
%     if(user.rateBased)% rate based selection for each segment
%         [val, index] = max(C_time < paramsDASH.l);
%         if(user.firstFlag)
%             user.rateBasedQual = paramsDASH.l(index - 1);
%             user.firstFlag = 0; 
%         end
%         user.currentQuality = user.rateBasedQual;
%     end
    %(buffer based decision) ENDS HERE    
    if (length(bThres) > 1) % DISCRETE quality levels
        if(abs(user.currentQuality - user.prevQuality) >= epsilon)
            user.qualitySwitches = user.qualitySwitches + 1;
        end
    else                    % CTX quality levels
        if(abs(user.currentQuality - user.prevQuality) >= paramsDASH.l(2) - paramsDASH.l(1))
            user.qualitySwitches = user.qualitySwitches + 1;
        end            
    end
%     user.squaredQualDiff = user.squaredQualDiff + (user.currentQuality - averageQualitySim) ^ 2 %computed from previous run
    user.avgQuality = user.avgQuality + user.currentQuality;%paramsDASH.l is in bits per slot!
    user.prevQuality = user.currentQuality;
    user.state = 1;            
end

if (user.state == 1)
    downProp = user.propOfSegmentDownloaded;
    remProp = 1 - downProp;
    %secsPerSeg = (1.0 / paramsDASH.segPerSec);
    slotPerSeg = (1.0 / segPerSlot);
    propSlot = -1;        
    %bitsRemainingforSegment = remProp * secsPerSeg * user.currentQuality;
    bitsRemainingforSegment = remProp * slotPerSeg * user.currentQuality;% current quality in bits per slot
    if (bitsRemainingforSegment < C_time) % C_time is the number of bits available in this simulation slot
        fprintf('%fbitsRemaining , %fC_time \n',bitsRemainingforSegment, C_time);
        downloaded = remProp;%number of segments downloaded in this slot            
        user.state = 0; %there are still some resources, DASH requests the next quality
        propSlot = bitsRemainingforSegment / C_time;
        user.slotsUsed = user.slotsUsed + propSlot;% rhs is the proportion of the slot used.
    else
        propSlot =  1;
        user.slotsUsed = user.slotsUsed + propSlot;
        user.propOfSegmentDownloaded = downProp + C_time / (user.currentQuality * slotPerSeg);
        downloaded = C_time / (user.currentQuality * slotPerSeg);%number of segments dowmloaded in this slot, it is mostly a fraction
    end
    if(user.inPrefetchStart == 1)% occurs at the very beginning of prefetching
        user.stateAtPrefetchStart = ScalarSystemState;%needs to be another field in the functions input
        user.inPrefetchStart = 0;
    end
    if (~user.isPlaying && (user.bufferLevel >= paramsDASH.qs))% the prefetching threshold is in terms of number of segments
        user.isPlaying = 1;
        % Find prefetch time:
        % For initial prefetching, starvedSegment is set to -1
        if (user.starvedSegment ==  -1)% -1 occurs only once, the prefetching time, that is what we initialized with
            user.prefetchTime = user.slotsUsed;% prefetch time in terms of number of slots used
            user.updateDelayMatrixFlag = 1;
            user.prefFlag = 0;
%             user.rateBased = 0;
        end;
    end
end

%printf('segments dowloaded in this slot = %d \n', downloaded);
% disp(downloaded);
played = segPerSlot * user.isPlaying;%we put 1 for the 2nd term, keeps playing for the whole slot
user.bufferLevel = user.bufferLevel + downloaded - played; %stores the bufferLevel (in number of segments)

if user.bufferLevel <= 0
    %starvation ; need to stop playing and prefetch qs segments
    user.bufferLevel = 0;
    user.isPlaying = 0;
    user.starvedSegment = user.currentSegment;
end
    
%else
% user.state = 2;
%update = 0;
end %function end
