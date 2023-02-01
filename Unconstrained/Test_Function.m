function [lb,ub,dim,fobj] = Test_Function(F)
dim = 30; % unless specified in the problem, dimension is set as 30 

switch F
    %%%%%%% %%%%%%%%%%%%unimodal
    
    %Sphere
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        
    %Powell Sum
    case 'F2'
        fobj = @F2;
        lb=-1;
        ub=1;
        
    %Schwefel 2.20
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        
    %Schwefel 2.21
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        
    %Step
    case 'F5'
        fobj = @F5;
        lb=-100;
        ub=100;
        
    %Schwefel 2.22
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        
    %Schwefel 2.23
    case 'F7'
        fobj = @F7;
        lb=-10;
        ub=10;
        
    %Rosenbrock
    case 'F8'
        fobj = @F8;
        lb=-30;
        ub=30;
        
    %Brown
    case 'F9'
        fobj = @F9;
        lb=-1;
        ub=4;
        
    %Dixon and Price
    case 'F10'
        fobj = @F10;
        lb=-10;
        ub=10;
        
    %Powell Singular
    case 'F11'
        fobj = @F11;
        lb=-4;
        ub=5;
        
    %Perm 0,D,Beta
    case 'F12'
        fobj = @F12;
        lb=-dim;
        ub=dim;
        
    %Sum Squares
    case 'F13'
        fobj = @F13;
        lb=-10;
        ub=10;
        
    %Schwefel 2.26
    case 'F14'
        fobj = @F14;
        lb=-500;
        ub=500;
        
    %Rastrigin
    case 'F15'
        fobj = @F15;
        lb=-5.12;
        ub=5.12;
        
    %Periodic
    case 'F16'
        fobj = @F16;
        lb=-10;
        ub=10;
        
    %Qing
    case 'F17'
        fobj = @F17;
        lb=-500;
        ub=500;
        
    %Alpine N. 1
    case 'F18'
        fobj = @F18;
        lb=-10;
        ub=10;
        
    %Xin-She Yang
    case 'F19'
        fobj = @F19;
        lb=-5;
        ub=5;
        
    %Ackley
    case 'F20'
        fobj = @F20;
        lb=-32;
        ub=32;
        
    %Trignometric 2
    case 'F21'
        fobj = @F21;
        lb=-500;
        ub=500;
        
    %Salomon
    case 'F22'
        fobj = @F22;
        lb=-100;
        ub=100;
        
    %Styblinski-Tang
    case 'F23'
        fobj = @F23;
        lb=-5;
        ub=5;
        
    %Griewank
    case 'F24'
        fobj = @F24;
        lb=-100;
        ub=100;
        
    %Xin-She Yang N. 4
    case 'F25'
        fobj = @F25;
        lb=-10;
        ub=10;
        
    %Xin-She Yang N. 2
    case 'F26'
        fobj = @F26;
        lb=-2*pi;
        ub=2*pi;
        
    %Gen. Pendlized
    case 'F27'
        fobj = @F27;
        lb=-50;
        ub=50;
        
    %Pendlized
    case 'F28'
        fobj = @F28;
        lb=-50;
        ub=50; 
        
    %Michalewics
    case 'F29'
        fobj = @F29;
        lb=0;
        ub=pi;    
            
    %Quartic Noise
    case 'F30'
        fobj = @F30;
        lb=-1.28;
        ub=1.28;
        
        %Shifted & Rotated Weierstrass     
        case 'F31' 
        fobj = @F31;
        lb=-100;
        ub=100;
        dim = 10;
        
        %Shifted & Rotated Happy Cat
        case 'F32'
        fobj = @F32;
        lb=-100;
        ub=100;
        dim = 10;
end
end

% F1 - Sphere
function o = F1(x)
    o=sum(x.^2);
end

% F2 - Powell Sum
function o = F2(x)
    n = size(x, 2);
    absx = abs(x);

    o = 0;
    for i = 1:n
        o = o + (absx(:, i) .^ (i + 1));
    end
end

% F3 - Schwefel 2.20
function o = F3(x)
	o = sum(abs(x), 2);
end

% F4 - Schwefel 2.21
function o = F4(x)
    o = max(abs(x), [], 2);
end

% F5 - Step
function o = F5(x)
    o=sum(abs((x+.5)).^2);
end

% F6 - Schwefel's 2.22
function o = F6(x)
    o=sum(abs(x))+prod(abs(x));
end

% F7 - Schwefel 2.23
function o = F7(x)
    o=sum(x .^10, 2);
end

% F8 - Rosenbrock
function o = F8(x)
    dim=size(x,2);
    o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F9 - Brown
function o = F9(x)
    dim=size(x,2);
    o = 0;
    for i = 1:dim-1
        xi1 = x(i+1)^2+1;
        xi = x(i)^2+1;
        o = o + (x(i)^2)^xi1 + (x(i+1)^2)^xi;
    end
end

% F10 - Dixon & Price
function o = F10(x)
    x1 = x(1);
    d = length(x);
    term1 = (x1-1)^2;
    sum = 0;
    for ii = 2:d
        xi = x(ii);
        xold = x(ii-1);
        new = ii * (2*xi^2 - xold)^2;
        sum = sum + new;
    end
    o = term1 + sum;
end

% F11 - Powell Singular
function o = F11(x)
    d = length(x);
    sum = 0;
    for ii = 1:(d/4)
        term1 = (x(4*ii-3) + 10*x(4*ii-2))^2;
        term2 = 5 * (x(4*ii-1) - x(4*ii))^2;
        term3 = (x(4*ii-2) - 2*x(4*ii-1))^4;
        term4 = 10 * (x(4*ii-3) - x(4*ii))^4;
        sum = sum + term1 + term2 + term3 + term4;
    end
    o = sum;
end

% F12 - Perm 0,D,Beta
function o = F12(x)
    if (nargin == 1)
        b = 10;
    end
    d = length(x);
    outer = 0;
    for ii = 1:d
        inner = 0;
        for jj = 1:d
            xj = x(jj);
    %         inner = inner + (jj^ii+b)*((xj/jj)^ii-1);
            inner = inner + (jj+b)*(xj^ii-(1/jj)^ii);
        end
        outer = outer + inner^2;
    end
    o = outer;
end

% F13 - Sum Squares
function o = F13(x)
    n = length(x);
    s = 0;
    for j = 1:n  
        s=s+j*x(j)^2; 
    end
    o = s;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%multilodal

% F14 - Schwefel's 2.26
function o = F14(x)
    dim=size(x,2);
    o=418.9829 - sum(-x.*sin(sqrt(abs(x))))/dim;
end

% F15 - Rastrigin
function o = F15(x)
    dim=size(x,2);
    o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F16 - Periodic
function o = F16(x)
    sin2x = sin(x) .^ 2;
    sumx2 = sum(x .^2, 2);
    o = 1 + sum(sin2x, 2) -0.1 * exp(-sumx2);
end

% F17 - Qing
function o = F17(x)
    n = size(x, 2);
    x2 = x .^2;
    
    o = 0;
    for i = 1:n
        o = o + (x2(:, i) - i) .^ 2;
    end
end

% F18 - Alpine N. 1
function o = F18(x)
     o = sum(abs(x .* sin(x) + 0.1 * x), 2);
end

% F19 - Xin-She Yang
function o = F19(x)
    n = size(x, 2);
    o = 0;
    for i = 1:n
        o = o + rand * (abs(x(:, i)) .^ i);
    end
end

% F20 - Ackley
function o = F20(x)
    n = size(x, 2);
    ninverse = 1 / n;
    sum1 = sum(x .^ 2, 2);
	sum2 = sum(cos(2 * pi * x), 2);
    
    o = 20 + exp(1) - (20 * exp(-0.2 * sqrt( ninverse * sum1))) - exp( ninverse * sum2);
end

% F21 - Trignometric 2
function o = F21(x)
    n = size(x, 2);
    
    o = 1;
    for i = 1:n
        o = o + 8*sin(7*(x(i)-0.9)^2)^2+6*sin(14*(x(1)-0.9)^2)^2+(x(i)-0.9)^2;
    end
end

% F22 - Salomon
function o = F22(x)
    x2 = x .^ 2;
    sumx2 = sum(x2, 2);
    sqrtsx2 = sqrt(sumx2);
    
    o = 1 - cos(2 .* pi .* sqrtsx2) + (0.1 * sqrtsx2);
end

% F23 - Styblinski-Tang
function o = F23(x)
    n = size(x, 2);
    o = 0;
    for i = 1:n
        o = o + ((x(:, i) .^4) - (16 * x(:, i) .^ 2) + (5 * x(:, i)));
    end
    o = 0.5 * o;
end

% F24 - Griewank
function o = F24(x)
    dim=size(x,2);
    o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F25 - Xin-She Yang N. 4
function o = F25(x)
     o = (sum(sin(x) .^2, 2) - exp(-sum(x .^ 2, 2))) .* exp(-sum(sin(sqrt(abs(x))) .^2, 2));
end

% F26 - Xin-She Yang N. 2
function o = F26(x)
    o = sum(abs(x), 2) .* exp(-sum(sin(x .^2), 2));
end

%Helping function for F27, F28
function o=Ufun(x,a,k,m)
    o=k.*((x-a).^m).*(x>a)+k.*((-x-a).^m).*(x<(-a));
end

% F27 - Gen. Pendlized
function o = F27(x)
    dim=size(x,2);
    o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
    ((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end

% F28 - Pendlized
function o = F28(x)
    dim=size(x,2);
    o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
    (1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F29 - Michalewics
function o = F29(x)
    n=size(x,2);
    m = 10;
    s = 0;
    for i = 1:n;
        s = s+sin(x(i))*(sin(i*x(i)^2/pi))^(2*m);
    end
    o = -s;
end

% F30 - Quartic Noise
function o = F30(x)
    dim=size(x,2);
    o=sum([1:dim].*(x.^4))+rand;
end

% F31 Shifted & Rotated Weierstrass
function o = F31(x)
    a =0.5;
    b =3.0;
    kMax = 20;
    dim = 10;
    sum = 0.0;
    shiftedMatrix = [4.4867071194977996e+01 8.6557399521842626e-01  -1.2297862364117918e+01   2.9827246270062048e+01   2.6528060932889602e+01  -6.2879900924339843e+01  -2.2494835379763892e+01   9.3017723082107295e+00   1.4887184097844738e+01  -3.1096867523666873e+01];
    rotatedMatrix = [
            -1.5433743057196678e-01 0.0000000000000000e+00 7.7666311726871273e-01 0.0000000000000000e+00 1.1571979400226866e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00;
             0.0000000000000000e+00 4.6806840267259536e-02 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00  -5.9264454599472804e-01 1.6314935476659614e-01 7.8737783169590370e-01 0.0000000000000000e+00 0.0000000000000000e+00;
            -1.7410812843826278e+00 0.0000000000000000e+00 -4.4194799352318298e-01 0.0000000000000000e+00 4.4605580480878959e-01 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00;
            0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 2.7077411154472419e-01 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00  -8.8999649318267127e-01 3.6686185770629254e-01;
            5.4888525059737507e-02 0.0000000000000000e+00 1.5570674387300532e+00 0.0000000000000000e+00  -3.0216546520289828e-01 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00;
            0.0000000000000000e+00 6.1164921138202333e-01 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00  -6.1748299284504526e-01  -1.5999277506278717e-01  -4.6797682388189477e-01 0.0000000000000000e+00 0.0000000000000000e+00;
            0.0000000000000000e+00 -1.1226733726002835e-01 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00  -1.3517591002752971e-01 9.4075663040175728e-01  -2.9000082877106131e-01 0.0000000000000000e+00 0.0000000000000000e+00;
            0.0000000000000000e+00 7.8172271740335475e-01 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 4.9921405128267116e-01 2.5052257846765580e-01 2.7736863877405393e-01 0.0000000000000000e+00 0.0000000000000000e+00;
            0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00  -3.0159372777109039e-01 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 2.8348126021733977e-01 9.1031840499614625e-01;
            0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 9.1417864987446651e-01 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 3.5713389257830902e-01 1.9165797370723034e-01
    ];
    %shifting
      shiftedX = x - shiftedMatrix;

    %rotating
    shitedRotaatedX = shiftedX * rotatedMatrix;

    for i =1: dim
        xi = shitedRotaatedX(1, i);
        inerSum_1 = 0.0;
        for k=0 :kMax
            inerSum_1 = inerSum_1 + (a^k * cos(2 * pi * b^k *(xi+0.5)));
        end
        sum = sum + inerSum_1;
    end
    inerSum_2 =0.0;
    for k=0 : kMax
        inerSum_2  = inerSum_2 + a^k * cos(pi * b^k);
    end
    sum = sum - dim * inerSum_2 +1;
    o = sum;
end

%F32 Shifted & Rotated Happy Cat
function o = F32(x)

  dim = 10;
  sum1=0.0;
  sum2=0.0;
  sum3=0.0;
  shiftedMatrix = [-6.0107960952496171e+00,  -6.3449972860258995e+01,  -3.6938623728667750e+00,  -2.7449007717635965e+00,  -5.3547271030744199e+01,   3.1015786282259867e+01,   2.3200459416583499e+00,  -4.6987858548289097e+01,   3.5061378905112562e+01,  -3.4047417731046465e+00];
  rotatedMatrix = [
                -7.6923624057192400e-02,  0.0000000000000000e+00,  7.2809258658661558e-02,  6.1371429917067155e-01,  0.0000000000000000e+00,  7.8239141541106805e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
                0.0000000000000000e+00, -1.1499983823069659e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  1.5729072158274271e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.3309066870600375e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
                -1.6730831752378217e-02,  0.0000000000000000e+00,  4.9480374519689890e-01,  6.5982384537901573e-01,  0.0000000000000000e+00, -5.6526261691115431e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
                9.1421044415115027e-01,  0.0000000000000000e+00, -3.3249140365486585e-01,  2.2489758522716782e-01,  0.0000000000000000e+00, -5.5586027556918202e-02,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
                0.0000000000000000e+00, -1.2704704488967578e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -7.6341623484218024e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  4.9980922801223232e-01,  0.0000000000000000e+00,  0.0000000000000000e+00;
                -3.9751993551989051e-01,  0.0000000000000000e+00, -7.9957334378299227e-01,  3.7068629354440513e-01,  0.0000000000000000e+00, -2.5544478964007222e-01,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
                0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  9.0184624725396623e-02,  0.0000000000000000e+00, -3.4243496198122719e-01, -9.3520320266563195e-01;
                0.0000000000000000e+00,  9.7696981452382070e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -6.8376531090322690e-01,  0.0000000000000000e+00,  0.0000000000000000e+00, -4.7094671586086240e-01,  0.0000000000000000e+00,  0.0000000000000000e+00;
                0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -6.1661688294908079e-01,  0.0000000000000000e+00, -7.5660067900409822e-01,  2.1757534831110126e-01;
                0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, -7.8208078427058847e-01,  0.0000000000000000e+00,  5.5704013261474550e-01, -2.7938492717261593e-01;
        ];
        %shifting
        ShiftedX = x - shiftedMatrix;
       
        %rotating
        shitedRotaatedX = ShiftedX * rotatedMatrix;

        for i= 1 :dim
            xi = shitedRotaatedX(1, i);
            sum1 = sum1  + xi^2- dim;
            sum2 = sum2 + xi^2;
            sum3 = sum3 + xi;
        end

        o =  (((abs(sum1))^1/4 )+(0.5 * sum2+sum3)/dim + 0.5);
end