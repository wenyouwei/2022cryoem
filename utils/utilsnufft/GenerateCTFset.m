function ctfs = GenerateCTFset(n,defocus)

if nargin == 1, defocus = [1.4 1.75 2];end
    
%     defocus=[2.035 1.4 1.9 1.86 1.875 2 1.58 1.45 1.75];
% defocus=[1.9 1.86 1.875 2 2.035 1.58 1.4 1.45 1.75];
alpha=0.07;res=3.36;lambda=0.0251;Cs=2;B=100;
% if mod(n,2)==0, siz=n+1;else siz = n; end
siz = n+1;
n_d  = length(defocus);  
ctfs = zeros(n,n,n_d);

for k=1:n_d
    %ctf = CTF(siz, res, lambda, defocus(k), Cs, B, alpha);
    ctf = CTF(siz, res, lambda, defocus(k), Cs, B, alpha);
       
    ctfs(:,:,k)=ctf(1:end-1,1:end-1);
    %if mod(n,2)==0, 
    %    ctfs(:,:,k) = ctf(1:end-1,1:end-1);
    %else
    %    ctfs(:,:,k) = ctf;
    %end
end


function h = CTF(n, pixA, lambda, defocus, Cs, B, alpha, deltadef, theta)
% f = CTF(n, pixA, lambda, defocus, Cs, B, alpha, deltadef, theta)
% f = CTF(n, pixA, Pars);  % Pick up the CTF parameters from the struct Pars
% f = CTF(n, Pars);  % use the Pars.pixA field too.
% Compute the contrast transfer function corresponding an n x n image with
% the sampling interval pixA in A/pixel.  Lambda is in A, Cs is in mm,
% B in A^2 and alpha in radians.  Defocus is in microns.  The last two
% arguments are optional.
% Alternatively, a struct Pars containing each of these fields is passed:
% pixA (optional); lambda, defocus, Cs, B, alpha, deltadef, theta.
%
% The CTF is computed with astigmatism if the struct (or argument list) has
% two additional fields, deltadef and theta.
%
% The result is returned in an nxn matrix with h(n/2+1,n/2+1) giving
% the zero-frequency amplitude.  Thus you must use ifftshift() on
% the result before performing a Fourier transform.  For example,
% to simulate the effect of the CTF on an image m, do this:
% fm=fftshift(fft2(m));
% cm=ifft2(ifftshift(fm.*ctf()));
%
% The first zero occurs at lambda*defocus*f0^2=1.
% e.g. when lambda=.025A, defocus=1um, then f0=1/16?A.

% Cs term fixed.  fs 4 Apr 04
% Struct option added fs 6 Jul 09
% Astig option added fs 18 Aug 09
%
astig=0;
cuton=0;
if isstruct(pixA)  % 2nd argument is a structure
    lambda=pixA;
    pixA=lambda.pixA;
end;
if isstruct(lambda) % 3rd argument is a structure
    P=lambda;
    lambda=P.lambda;
    defocus=P.defocus;
    Cs=P.Cs;
    B=P.B;
    alpha=P.alpha;
    if isfield(P,'deltadef')  % we are handling astigmatism
        deltadef=P.deltadef;
        theta=P.theta;
        if deltadef ~= 0
            astig=1;
        end;
    end;
        if isfield(P,'cuton')
        cuton=P.cuton;
    end;

elseif nargin>7 % we have astig parameters
        astig=1;
end;

if astig
    [r1 ang]=RadiusNorm(n,fctr(n));
    % we use the defocus formula from Henderson '86:
    df=defocus+deltadef*cos(2*(ang-theta));
else
    r1=RadiusNorm(n,fctr(n));
    df=defocus;
end;
r2=r1.^2;
f0 = 1./pixA;  % Spatial frequency unit (inverse ?)

k2=-df*pi*lambda*1e4*f0.^2;  % this may be a matrix
k4= pi/2*Cs*lambda^3*1e7*f0.^4;  % this is a scalar.
kr=f0^2*B;  % B-factor

if Cs==0
    h=sin(k2.*r2-alpha).*exp(-kr*r2);
else
    h=sin(k2.*r2+k4*r2.*r2-alpha).*exp(-kr*r2);
end

if cuton  % handle sharp cut-on of a phase plate.
    h=h.*(0.5+0.5*erf((r1*f0-cuton)*10/cuton));
end;

function org=fctr(n)
% function org=fctr(n)
% Center of an FFT-shifted image.  We use this center coordinate for all
% rotations and centering operations.
% If n is even, org=[n/2+1;n/2+1].
% If n is odd, org is [(n+1)/2;(n+1)/2].
% n can be a two-element vector n=[nx ny]
if numel(n)<2
    n=[n n];
end;
org=ceil((n+1)/2);


function [r theta]=RadiusNorm(n,org)
% function [r theta]=RadiusNorm(n,org)
% Create an n(1) x n(2) array where the value at (x,y) is the distance from the
% origin, normalized such that a distance equal to the width or height of
% the array = 1.  This is the appropriate function to define frequencies
% for the fft of a rectangular image.  
% For a square array of size n (or [n n]) the following is true:
% RadiusNorm(n) = Radius(n)/n.
% If a second output argument is given,
% the angle (in radians) is also supplied.
% The org argument is optional, in which case the FFT center is used.
% 
if numel(n)<2
    n=[n n];
end;
if nargin<2
     org=ceil((n+1)/2);
end;

n=single(n);  % this will force the output to be single
x=org(1);
y=org(2);
[X,Y]=ndgrid((1-x:n(1)-x)/n(1),(1-y:n(2)-y)/n(2));  % Make zero at x,y
    r=sqrt(X.^2+Y.^2);
if nargout>1
    theta=atan2(Y,X);
end;

