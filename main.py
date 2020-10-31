#-------------------------------------------------------------------------------
#
# Autores : Isaac Macedo, Gustavo Nogueira
# Data: Setembro/2020
# 
# Licença:
# Este programa é um software livre; você pode redistribuí-lo e/ou
# modificá-lo sob os termos da Licença Pública Geral GNU como publicada
# pela Free Software Foundation; na versão 3 da Licença, ou
# (a seu critério) qualquer versão posterior.
#
# Este programa é distribuído na esperança de que possa ser útil,
# mas SEM NENHUMA GARANTIA; sem uma garantia implícita de ADEQUAÇÃO
# a qualquer MERCADO ou APLICAÇÃO EM PARTICULAR. Veja a
# Licença Pública Geral GNU para mais detalhes.
#
# Você deve ter recebido uma cópia da Licença Pública Geral GNU junto
# com este programa. Se não, veja <http://www.gnu.org/licenses/>.
#
# Uso:
# 	python3 main.py -i arquivoDeEntrada, onde i = 1,2,3,4
#-------------------------------------------------------------------------------
import sys
import os.path
import subprocess
import distutils.spawn
#-------------------------------------------------------------------------------
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------------
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
#-------------------------------------------------------------------------------

def plotImage(name, data, levels, wm, wl, sigmax, sigmay, texto, nome):
	N = 1000 + 1j
	x,y,z=np.genfromtxt(data, unpack=True)

	grid_x, grid_y = np.mgrid[x.min():x.max():N, y.min():y.max():N]
	#grid_z=s.griddata((x,y), z, (grid_x, grid_y), method='cubic')
	#grid_z=s.griddata((x,y), z, (grid_x, grid_y), method='nearest')
	grid_z=griddata((x,y), z, (grid_x, grid_y), method='linear')

	plt.contour(grid_x, grid_y, grid_z, levels=levels, colors='#000000')
	#plt.contourf(grid_x, grid_y, grid_z, levels=levels, colors=['#055698', '#1E68A7', '#5D8DBD', '#FFFFFF'])

	# color ->
	# b: blue, g: green, r: red, c: cyan, m: magenta, y: yellow, k: black, w: white
	# marker ->
	# .: point, ,: pixel, o: circle, v: triangle_down, ^: triangle_up, >: triangle_right
	# <: triangle_left, 1: tri_down, 2: tri_up, 3: tri_left, 4: tri_right, s: square
	# p: pentagon, *: star, h: hexagon1, H: hexagon2, +: plus, x: x, D: diamon, 
	# d: thin_diamond, |, vline, _: hline  
	plt.plot(wm, wl, 'k*')
	plt.title(nome)
	plt.xlabel(r'$\Omega_m$', {'color': 'k', 'fontsize': 14})
	plt.ylabel(r'$\Omega_\Lambda$', {'color': 'k' ,'fontsize': 14})
	plt.text(.7, .9, r'$\Omega_m$: {}$\pm${:.3f}'.format(wm, sigmax), {'color': 'k', 'fontsize': 12})
	plt.text(.7, .85, r'$\Omega_\Lambda$: {}$\pm${:.3f}'.format(wl, sigmay), {'color': 'k', 'fontsize': 12})

	plt.text(.7, .1, texto, {'color': 'k', 'fontsize': 12})

	#plt.colorbar()
	#plt.grid(True)
	#plt.show()
	plt.savefig(name)
	
#-------------------------------------------------------------------------------

def textPlot(wmmin, wlmin, wkmin, Mumin, hmin, c6min, levels1, levels2, levels3, NDOF, param):

	t01=r"""
$\nu$ = {NDOF}
$\Omega_k$ = {wkmin:.3f}""".format(NDOF=NDOF,  
		wkmin=wkmin)
	
	t02=r"""
$\mu$ = {Mumin:.3f}
$h_m$ = {hmin:.3f}""".format(Mumin=Mumin, 
		hmin=hmin)
		
	t04=r"""
$\chi^2/\nu$ = {chi_nu:.3f}
1$\sigma$ = {prob68:.3f}
2$\sigma$ = {prob95:.3f}
3$\sigma$  = {prob99:.3f}""".format(chi_nu=c6min/NDOF, 
		prob68=levels1, 
		prob95=levels2,
		prob99=levels3)
		
	if param == "1":
		return t01+t02+t04
	elif param == "2":
		return t01+t04
	elif param == "3" or param == "4":
		return t01+t04
	elif param == "5":
		return t01+t02
		
#-------------------------------------------------------------------------------
		
def main():

	lista_arg = sys.argv
	txt = """
Uso: programa.py parametro arquivo_entrada"
Parametros: 
	-1: hdl
	-2: mag
	-3: da wang
	-4: riess
	-5: cmbBAO
"""
	
	if len(lista_arg) < 3:
		print(txt)
		exit(1)

	param = {"-1":"1", "-2":"2", "-3":"3", "-4":"4", "-5":"5"}.get(lista_arg[1], "emp")
	if param == "emp":
		print(txt)
		return

	compiler = "gfortran"
	if not distutils.spawn.find_executable(compiler):
		print("Compilador \'", compiler, "\' não encontrado!")
		exit(1)
			
	if param == "5":
		filesIn = [lista_arg[2], lista_arg[3]]
	else:
		filesIn = [lista_arg[2]]
		
	fileFortran = ["main_mod2.f95", "gamma.f95", "integral.f95", "derivada.f95"]+filesIn
	
	for file in fileFortran:
		if not os.path.isfile(file):
			print("Arquivo \'",file," \'não encontrado!")
			exit(1)
	
	file_exe = "file.exe"
	if not distutils.spawn.find_executable(file_exe):
		if param == "5":
			subprocess.run([compiler, "-o", file_exe]+fileFortran[:-2])
		else:
			subprocess.run([compiler, "-o", file_exe]+fileFortran[:-1])
	
	map_param = {"1":"hdl-", "2":"mag-", "3":"daWang-", "4":"riess-", "5":"cmbBAO-"}
	
	if param == "5":
		arquivoDeEntrada = [lista_arg[2], lista_arg[3]]
	else:
		arquivoDeEntrada = [lista_arg[2]]
		
	arquivoDeSaida = map_param.get(param)[:-1]+".dat"#+arquivoDeEntrada
	arquivoSaidaPython = map_param.get(param)+"Python.dat"#+arquivoDeEntrada
	
	run_proc = ["./"+file_exe, param] + arquivoDeEntrada + [arquivoDeSaida, arquivoSaidaPython]
	subprocess.run(run_proc)
	
	filePython = np.loadtxt(arquivoSaidaPython, usecols=1)

	NDOF, wmmin, wlmin, wkmin, Mumin, hmin, c6min, chi2, prob68, prob95, prob99, err1, err2 = filePython
		
	textoPlot = textPlot(wmmin, wlmin, wkmin, Mumin, hmin, c6min, prob68, prob95, prob99, NDOF, param)
						
	plotImage(arquivoDeSaida[:-4], arquivoDeSaida,
			 [0, prob68, prob95, prob99], wmmin, wlmin, err1, err2, textoPlot, map_param.get(param))
			
	return
			
if __name__ == "__main__":
	main()
