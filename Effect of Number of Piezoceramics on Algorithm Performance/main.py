import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import os 


def main():
    os.system("cls")
    
    excel_filename = "n_10 estimation results.xlsx"
    
    df = pd.read_excel(excel_filename)

    # all values in mL
    sphere_means = np.asarray([df.iloc[22, 1], df.iloc[22, 5], df.iloc[22, 9], df.iloc[22, 13], df.iloc[51, 1], df.iloc[51, 5], df.iloc[51, 9]])
    sphere_errors = np.asarray([df.iloc[24, 1], df.iloc[24, 5], df.iloc[24, 9], df.iloc[24, 13], df.iloc[53, 1], df.iloc[53, 5], df.iloc[53, 9]])
    
    ellipse_means = [df.iloc[22, 2], df.iloc[22, 6], df.iloc[22, 10], df.iloc[22, 14], df.iloc[51, 2], df.iloc[51, 6], df.iloc[51, 10]]
    ellipse_errors = [df.iloc[24, 2], df.iloc[24, 6], df.iloc[24, 10], df.iloc[24, 14], df.iloc[53, 2], df.iloc[53, 6], df.iloc[53, 10]]
    
    x_axis = [4,5,6,7,8,9,10]
    
    plt.figure()
    # plt.title("Effect of number of coordinates on spherical fitting algorithm performance\n(spherical-shaped round-bottom flask)")
    # plt.title("Spherical-shaped round-bottom flask")
    plt.errorbar(x_axis, sphere_means, yerr=sphere_errors, fmt=":o", lw=1.5, color="red", ecolor="black", capsize=10, 
                    capthick=2, barsabove=True)
    plt.xlabel("Number of Piezoceramic Elements"); plt.ylabel("Estimated Volume [mL]")
    
    plt.savefig("spherical flask - effect of coord number with error bars.PNG", dpi=300, bbox_inches="tight")
    plt.savefig("spherical flask - effect of coord number with error bars.SVG", dpi=300, bbox_inches="tight")
    
    




    plt.figure()
    # plt.title("Effect of number of coordinates on spherical fitting algorithm performance\n(oval-shaped round-bottom flask)")
    # plt.title("Oval-shaped round-bottom flask")
    plt.errorbar(x_axis, ellipse_means, yerr=ellipse_errors, fmt=":o", lw=1.5, color="red", ecolor="black", capsize=10, 
                    capthick=2, barsabove=True)
    plt.xlabel("Number of Piezoceramic Elements"); plt.ylabel("Estimated Volume [mL]")
    
    plt.savefig("oval flask - effect of coord number with error bars.PNG", dpi=300, bbox_inches="tight")
    plt.savefig("oval flask - effect of coord number with error bars.SVG", dpi=300, bbox_inches="tight")
    


if __name__ == "__main__":
    main()
    plt.show()