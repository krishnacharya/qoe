function [user,n_users] = DASH_SP(paramsDASH, segPerSlot, user, time, C_time, n_users, bThres, idThres)
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
totalC = C_time;
epsilon = 1e-5;
propSlot = 0;
if user.currentSegment <= user.totalNbrOfSegments
    while (C_time > 0)
        if (user.state == 0)
            if user.currentSegment == user.totalNbrOfSegments
                user.state = 2;
                user.avgQuality = user.avgQuality/user.totalNbrOfSegments;
                user.exitTime = time;
                if (user.prefetchTime <= 0) % Video was delivered completely in prefetching.
                    user.prefetchTime = time;
                end
                n_users = n_users-1;
                return;
            end
            user.currentSegment = user.currentSegment + 1;
            
            % Make decision on the next quality
            % (buffer based decision)  STARTS HERE
            if user.bufferLevel <= paramsDASH.Bmin
                user.currentQuality=paramsDASH.l(1);
            else
                if (length(bThres) > 1) % DISCRETE quality levels
                    idx = max(idThres(bThres <= user.bufferLevel));
                    user.currentQuality = paramsDASH.l(idx);
                else                    % CTX quality levels
                    user.currentQuality = paramsDASH.rate*user.bufferLevel + paramsDASH.c;
                    user.currentQuality = min(paramsDASH.l(end), user.currentQuality);
                end
            end
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
            user.avgQuality = user.avgQuality + user.currentQuality;
            user.prevQuality = user.currentQuality;
            user.state = 1;
            
        end
        
        if (user.state == 1)
            prev = user.propOfDownloadedSegment;
            capacityUsed = user.currentQuality/segPerSlot * (1 - prev);
            if (capacityUsed < C_time)
                C_time = C_time - capacityUsed;
                downloaded = (1 - prev);
                user.propOfDownloadedSegment = 0;
                user.state = 0; %there are still some resources, DASH requests the next quality
                propSlot = capacityUsed/totalC;
            else
                propSlot=C_time/totalC;
                user.propOfDownloadedSegment = prev +  ...
                    C_time * segPerSlot / user.currentQuality;
                downloaded = user.propOfDownloadedSegment - prev;
                C_time = 0;
            end
            if (~user.isPlaying && user.bufferLevel>= paramsDASH.qs)
                user.isPlaying=1;
                % Find prefetch time:
                % For initial prefetching, starvedSegment is set to -1
                if (user.starvedSegment ==  -1)
                    user.prefetchTime = time + propSlot;
                end;
            end
        end
        
        played = segPerSlot * propSlot * user.isPlaying ;
        user.bufferLevel = user.bufferLevel + downloaded - played ;
        if user.bufferLevel <= 0
            %starvation ; need to stop playing and prefetch qs segments
            user.bufferLevel = 0;
            user.isPlaying = 0;
            user.starv = 1;
            user.starvedSegment = user.currentSegment;
        end
    end
    
else
    user.state = 2;
    user.exitTime = time;
end
end
