/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * @author: Gabriel
 * @date:   27.04.2017
 */

namespace anpi {

  template <typename T>
  PlotTNSHA<T>::PlotTNSHA(){}

  template <typename T>
  void PlotTNSHA<T>::initialize(){
    Py_Initialize();
    PyRun_SimpleString("import matplotlib.pyplot as plt");
    PyRun_SimpleString("import numpy as np");
  }

  template <typename T>
  void  PlotTNSHA<T>::plot(std::vector<T>& datax,std::vector<T>& datay,std::vector<T>& datau,std::vector<T>& datav) {


    // Convert the vectors of data into Python strings
    std::string xstr  = "X = [";
    std::string ystr  = "Y = [";
    std::string ustr  = "U = [";
    std::string vstr  = "V = [";

    char c=',';
    for(size_t i = 0; i < datax.size(); i++) {
          if (i == datax.size()-1) {
            c=']';
          }
    xstr.append(std::to_string(datax[i])   + c);
    ystr.append(std::to_string(datay[i])   + c);

  }
   c=',';
   for(size_t i = 0; i < datav.size(); i++) {
            if (i == datav.size()-1) {
              c=']';
            }
      ustr.append(std::to_string(datau[i])   + c);
      vstr.append(std::to_string(datav[i])   + c);

   }


    PyRun_SimpleString(xstr.c_str());
    PyRun_SimpleString(ystr.c_str());
    PyRun_SimpleString(ustr.c_str());
    PyRun_SimpleString(vstr.c_str());

    PyRun_SimpleString("fig, ax = plt.subplots()");
    PyRun_SimpleString("M = np.hypot(U, V)");
    PyRun_SimpleString("ax.quiver(X, Y, U, V,M)");
    PyRun_SimpleString("plt.plot(X, Y, linewidth=5, color='black')");


}

  template <typename T>
  void PlotTNSHA<T>::show(){
    PyRun_SimpleString("plt.show()");
  }

} // namespace anpi