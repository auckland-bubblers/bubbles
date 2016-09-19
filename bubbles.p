# GnuPlot script file bubbles.p
#
# gnuplot> cd 'storage/sdcard0/Documents/Research/Gnuplot/Bubbles'
# gnuplot> load "bubbles.p"
#
unset multiplot
set key noautotitle
#
#
# Input parameters
#
alpha_phi = 0.1
m_phi = 1
phi_0 = 1
tau_star = 10
x_p = 128
n = 1
#
#
# Fresnel integrals
#
S(x) = x < 2.75 ? ( sum [i=0:11] (-1)**i * x**(4*i+3) / ((4*i+3)*(2*i+1)!) ) : ( sqrt(pi/8) - cos(x**2)/(2*x) - sin(x**2)/(4*x**3) + 3*cos(x**2)/(8*x**5) + 3*5*sin(x**2)/(16*x**7) )
#
C(x) = x < 2.75 ? ( sum [i=0:11] (-1)**i * x**(4*i+1) / ((4*i+1)*(2*i)!) ) : ( sqrt(pi/8) + sin(x**2)/(2*x) - cos(x**2)/(4*x**3) - 3*sin(x**2)/(8*x**5) + 3*5*cos(x**2)/(16*x**7) )
#
set title "Fresnel integrals"
#
set output "S.png"
set xlabel "x"
set ylabel "S(x)"
set xrange [0:10]
set yrange [0:1]
set samples 1e4
plot S(x) lc rgb "#0000FF"
#
set output "C.png"
set xlabel "x"
set ylabel "C(x)"
set xrange [0:10]
set yrange [0:1]
set samples 1e4
plot C(x) lc rgb "#0000FF"
#
#
# Analytic bubble solution
#
phi(tau) = tau < n/(2*m_phi) ? 0 : ( tau < tau_star ? ( (alpha_phi*phi_0/(1+alpha_phi)) * (tau_star/tau)**(n/2) * exp(m_phi*(tau-tau_star)) ) : ( phi_0 - (phi_0/sqrt(1+alpha_phi)) * sqrt(tau_star/tau)**n * cos(sqrt(alpha_phi)*m_phi*(tau-tau_star)+atan(sqrt(alpha_phi))) ) )
#
set terminal png
set output "phi.png"
set tmargin 2
set title ""
set multiplot title "Analytic bubble solution for: alpha_{phi} = " . gprintf("%g",alpha_phi) . " , tau_* = " . gprintf("%g",tau_star)
set xlabel "x"
set ylabel "phi(x)"
set xrange [0:x_p]
set yrange [0:2*phi_0]
set samples 1e4
t = tau_star
plot phi(sqrt(t**2-x**2)) lc rgb "#FF0000"
t = sqrt((1*x_p/4)**2+tau_star**2)
plot phi(sqrt(t**2-x**2)) lc rgb "#FFFF00"
t = sqrt((2*x_p/4)**2+tau_star**2)
plot phi(sqrt(t**2-x**2)) lc rgb "#00FF00"
t = sqrt((3*x_p/4)**2+tau_star**2)
plot phi(sqrt(t**2-x**2)) lc rgb "#00FFFF"
t = sqrt(x_p**2+tau_star**2)
plot phi(sqrt(t**2-x**2)) + ( x < 2*x_p-t ? 0 : phi(sqrt(t**2-(2*x_p-x)**2)) ) lc rgb "#0000FF"
unset tmargin
unset multiplot
#
#
# Analytic spectrum parameters and functions for n=1
#
t = sqrt(x_p**2+tau_star**2)
#
a = atan(sqrt(alpha_phi))
b = asin(tau_star/t)
m = sqrt(alpha_phi)*m_phi
#
k(l) = pi*l/x_p
w(l) = sqrt(k(l)**2+m**2)
g(l) = atan(m/k(l))
#
eta(l) = atan(k(l)*tan(b)/m_phi)
#
f_ro(l) = 2*(phi_0/(k(l)*x_p)) * sin(a)**2 * sin(eta(l)) * cos(k(l)*t*cos(b)+eta(l))
#
fdot_ro(l) = -2*k(l)*(phi_0/(k(l)*x_p)) * sin(a)**2 * (sin(eta(l))/cos(b)) * ( sin(k(l)*t*cos(b)+eta(l)) + (cos(eta(l))/(k(l)*t*cos(b)))*cos(k(l)*t*cos(b)+2*eta(l)) )
#
f_ph(l) = 2*(phi_0/(k(l)*x_p)) * sin(k(l)*t*cos(b))
#
fdot_ph(l) = 2*k(l)*(phi_0/(k(l)*x_p)) * (cos(k(l)*t*cos(b))/cos(b))
#
Asc(l) = w(l)*t*(sin(b-g(l)))**2/(2*cos(b-g(l)))
#
A_os(l) = g(l) > b ? ( sqrt((sin(b))**(n-1)*(sin(g(l)))**(2-n)) * ( cos(w(l)*t-m*tau_star+a) * ( sqrt(pi/8) + C(sqrt(w(l)*t/2)*(g(l)-b)) ) + sin(w(l)*t-m*tau_star+a) * ( sqrt(pi/8) + S(sqrt(w(l)*t/2)*(g(l)-b)) ) ) ) : ( sqrt(sin(b)/cos(b-g(l))) * ( cos(Asc(l)+k(l)*t*cos(b)+a) * ( sqrt(pi/8) - C(sqrt(Asc(l))) ) + sin(Asc(l)+k(l)*t*cos(b)+a) * ( sqrt(pi/8) - S(sqrt(Asc(l))) ) ) )
#
Adot_os(l) = g(l) > b ? ( sqrt((sin(b))**(n-1)*(sin(g(l)))**(2-n)) * ( -w(l)*sin(w(l)*t-m*tau_star+a) * ( sqrt(pi/8) + C(sqrt(w(l)*t/2)*(g(l)-b)) ) + w(l)*cos(w(l)*t-m*tau_star+a) * ( sqrt(pi/8) + S(sqrt(w(l)*t/2)*(g(l)-b)) ) + sqrt(w(l)/(2*t)) * (tan(b)+(g(l)-b)/2) * cos(w(l)*t-m*tau_star+a-(g(l)-b)**2*w(l)*t/2) ) ) : ( (w(l)/(2*cos(b))) * sqrt(sin(b)/(cos(b-g(l)))**5) * (2*cos(b)*cos(b-g(l))-cos(g(l))*(sin(b-g(l)))**2) * ( -sin(Asc(l)+k(l)*t*cos(b)+a) * ( sqrt(pi/8) - C(sqrt(Asc(l))) ) + cos(Asc(l)+k(l)*t*cos(b)+a) * ( sqrt(pi/8) - S(sqrt(Asc(l))) ) ) + (1/(2*cos(b))) * sqrt(w(l)/(2*t)) * sqrt(sin(b))/cos(b-g(l))**2 * (sin(b)+sin(g(l))*cos(b-g(l))) * cos(k(l)*t*cos(b)+a) )
#
Bsc(l) = w(l)*t*(sin(b+g(l)))**2/(2*cos(b+g(l)))
#
B_os(l) = sqrt(sin(b)/cos(b+g(l))) * ( cos(Bsc(l)+k(l)*t*cos(b)-a) * ( sqrt(pi/8) - C(sqrt(Bsc(l))) ) + sin(Bsc(l)+k(l)*t*cos(b)-a) * ( sqrt(pi/8) - S(sqrt(Bsc(l))) ) )
#
Bdot_os(l) = (w(l)/(2*cos(b))) * sqrt(sin(b)/(cos(b+g(l)))**5) * (2*cos(b)*cos(b+g(l))-cos(g(l))*(sin(b+g(l)))**2) * ( -sin(Bsc(l)+k(l)*t*cos(b)-a) * ( sqrt(pi/8) - C(sqrt(Bsc(l))) ) + cos(Bsc(l)+k(l)*t*cos(b)-a) * ( sqrt(pi/8) - S(sqrt(Bsc(l))) ) ) + (1/(2*cos(b))) * sqrt(w(l)/(2*t)) * sqrt(sin(b))/cos(b+g(l))**2 * (sin(b)-sin(g(l))*cos(b+g(l))) * cos(k(l)*t*cos(b)-a)
#
f_os(l) = -k(l)*(phi_0/(k(l)*x_p))*sqrt(2*tau_star/w(l))*cos(a)*( A_os(l) + B_os(l) )
#
fdot_os(l) = -k(l)*(phi_0/(k(l)*x_p))*sqrt(2*tau_star/w(l))*cos(a)*( Adot_os(l) + Bdot_os(l) )
#
f(l) = f_ro(l) + f_ph(l) + f_os(l)
#
fdot(l) = fdot_ro(l) + fdot_ph(l) + fdot_os(l)
#
#
# Analytic spectrum graphs
#
set title "Spectrum for: alpha_{phi} = " . gprintf("%g",alpha_phi) . " , tau_* = " . gprintf("%g",tau_star) . " , x_p = " . gprintf("%g",x_p)
#
set terminal png
set output "f_ro.png"
set xlabel "l"
set ylabel "k(l)*f_{ro}(l)"
set xrange [1:1e4]
set yrange [*:*]
set samples 1e4
set logscale x
plot k(x)*f_ro(x) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "fdot_ro.png"
set xlabel "l"
set ylabel "fdot_{ro}(l)"
set xrange [1:1e4]
set yrange [*:*]
set samples 1e4
set logscale x
plot fdot_ro(x) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "f_ph.png"
set xlabel "l"
set ylabel "k(l)*f_{ph}(l)"
set xrange [1:1e1]
set yrange [*:*]
set samples 1e4
set logscale x
plot k(x)*f_ph(x) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "fdot_ph.png"
set xlabel "l"
set ylabel "fdot_{ph}(l)"
set xrange [1:1e1]
set yrange [*:*]
set samples 1e4
set logscale x
plot fdot_ph(x) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "f_os.png"
set xlabel "l"
set ylabel "k(l)*f_{os}(l)"
set xrange [1:1e5]
set yrange [*:*]
set samples 1e4
set logscale x
plot k(x)*f_os(x) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "fdot_os.png"
set xlabel "l"
set ylabel "fdot_{os}(l)"
set xrange [1:1e5]
set yrange [*:*]
set samples 1e4
set logscale x
plot fdot_os(x) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "f.png"
set xlabel "l"
set ylabel "k(l)*f(l)"
set xrange [1:1e4]
set yrange [*:*]
set samples 1e4
set logscale x
plot k(x)*f(x) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "fdot.png"
set xlabel "l"
set ylabel "fdot(l)"
set xrange [1:1e4]
set yrange [*:*]
set samples 1e4
set logscale x
plot fdot(x) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "ke.png"
set xlabel "l"
set ylabel "k(l)*fdot(l)^2"
set xrange [1:1e4]
set yrange [0:*]
set samples 1e4
set logscale x
plot k(x)*fdot(x)**2 lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "ge.png"
set xlabel "l"
set ylabel "k(l)^3*f(l)^2"
set xrange [1:1e4]
set yrange [0:*]
set samples 1e4
set logscale x
plot k(x)**3*f(x)**2 lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "pe.png"
set xlabel "l"
set ylabel "k(l)*mu^2*f(l)^2"
set xrange [1:1e4]
set yrange [0:*]
set samples 1e4
set logscale x
plot k(x)*m**2*f(x)**2 lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "e.png"
set xlabel "l"
set ylabel "k(l)*E(l)"
set xrange [1:1e4]
set yrange [0:*]
set samples 1e4
set logscale x
plot k(x)*(fdot(x)**2+(k(x)**2+m**2)*f(x)**2) lc rgb "#0000FF"
unset logscale
#
set terminal png
set output "f_ro,f_ph,f_os.png"
set tmargin 2
set title ""
set multiplot title "Spectrum for: alpha_{phi} = " . gprintf("%g",alpha_phi) . " , tau_* = " . gprintf("%g",tau_star) . " , x_p = " . gprintf("%g",x_p)
set xlabel "l"
set ylabel "k(l)*f(l)"
set xrange [1:1e4]
set yrange [*:*] writeback
set samples 1e3
set logscale x
plot k(x)*f_os(x) lc rgb "#0000FF"
set yrange restore
plot k(x)*f_ph(x) lc rgb "#00FF00"
plot k(x)*f_ro(x) lc rgb "#FF0000"
unset tmargin
unset multiplot
unset logscale
#
