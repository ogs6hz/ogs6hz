z = 10000; %z distance 
b = 150; %radius of length
vel = 1540;  %speed of sound
a = 100; %radius of width
xvalues = [140, 100, 60, 60, 0]; %x-direction of field point
yvalues = [187.5, 112.5, 112.5, 0, 0]; %y-direction of field point
xsize = length(xvalues);
tvalues = 6.4934:1/100e6:6.5;
tsize = length(tvalues);
figure
for i = 1:xsize
    x = abs(xvalues(i));
    y = abs(yvalues(i));
    h = [];
    if abs(x) <= a && abs(y) <= b
        small = [a+x,a-x,a+x,a-x];
        elmall = [b+y,b+y,b-y,b-y];
        coeffall = [1,1,1,1];
    elseif abs(x) > a && abs(y) > b
        small = [2*a,2*a,x-a,x-a];
        elmall = [2*b,y-b,2*b,y-b];
        coeffall = [1,-1,-1,-1];
    elseif abs(x) > a && abs(y) <= b 
        small = [2*a,2*a,x-a,x-a];
        elmall = [b+y,b-y,b+y,b-y];
        coeffall = [1,1,-1,-1];
    elseif abs(x) <= a && abs(y) > b 
        small = [a+x,a+x,a-x,a-x];
        elmall = [2*b,y-b,2*b,y-b];
        coeffall = [1,-1,1,-1];
    end
    for k = 1:tsize
        t = tvalues(k); 
        t1 = z/vel;
        summation = 0;
        for m = 1:4
            sm = small(m);
            elm = elmall(m);
            coeff = coeffall(m);
            t2 = (1/vel)*sqrt(z^2 + sm^2);
            t3 = (1/vel)*sqrt(z^2 + elm^2);
            t4 = (1/vel)*sqrt(z^2 + sm^2 + elm^2);
            if t<t1
                var = 0;
            elseif t<t2
                var = pi/2;
            elseif t<t3
                var = pi/2 - acos(sm./sqrt((vel.^2)*(t.^2)-(z.^2)));
            elseif t<=t4
                var = pi/2 - acos(sm./sqrt((vel.^2)*(t.^2)-(z.^2))) - acos(elm./sqrt((vel.^2).*(t.^2)-(z.^2)));
            else
                var = 0;
            end                
            if var ~= 0
                summation = summation + coeff*var; 
            end
        end
        h(k) = (vel/(2*pi))*summation;  
    end
    subplot(xsize,1,i)
    plot(tvalues,h)
    title(['Impulse Response at z/a = 10, x/a = ', num2str(x/a), ' and y/b = ', num2str(y/b)])
    xlabel('t')
    ylabel('h(x,t)')
    %axis tight
    xlim([6.493 6.5])
    %ylim([0 1540])
    hold on
end