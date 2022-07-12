/*
   This file is part of FERRET, an add-on module for MOOSE

   FERRET is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>
   and be sure to track new changes at github.com/mangerij/ferret

**/

#include "FerretApp.h"

//Moose Includes
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

// Create a performance log
PerfLog Moose::perf_log("Ferret");

 // Begin the main program.
int main(int argc, char *argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // Register this application's MooseApp and any it depends on
  FerretApp::registerApps();

  // This creates dynamic memory that we're responsible for deleting
  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("FerretApp", argc, argv);

  std::cout<<"   .-.                                                          ___       "<<"\n";
  std::cout<<"  /    \\                                                       (   )     "<<"\n";
  std::cout<<"  | .`. ;      .--.      ___           ___            .--.      | |_      "<<"\n";
  std::cout<<"  | |(___)    /    \\    (   )         (   )          /    \\    (   __)  "<<"\n";
  std::cout<<"  | |_       |  .-. ;    | ' .-. ;     | ' .-. ;    |  .-. ;    | |       "<<"\n";
  std::cout<<" (   __)     |  | | |    |  / (___)    |  / (___)   |  | | |    | | ___   "<<"\n";
  std::cout<<"  | |        |  |/  |    | |           | |          |  |/  |    | |(   )  "<<"\n";
  std::cout<<"  | |        |  ' _.'    | |           | |          |  ' _.'    | | | |   "<<"\n";
  std::cout<<"  | |        |  .'.-.    | |           | |          |  .'.-.    | ' | |   "<<"\n";
  std::cout<<"  | |        '  `-' /    | |           | |          '  `-' /    ' `-' ;   "<<"\n";
  std::cout<<" (___)        `.__.'    (___)         (___)          `.__.'      `.__.    "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"            ...a MOOSE package for simulating the ferroic nanostructure   "<<"\n";

  std::cout<<"                                                                          "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"__________________________________________________________________________"<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"   FERRET is free software: you can redistribute it and/or modify         "<<"\n";
  std::cout<<"   it under the terms of the GNU General Public License as published by   "<<"\n";
  std::cout<<"   the Free Software Foundation, either version 3 of the License, or      "<<"\n";
  std::cout<<"   (at your option) any later version.                                    "<<"\n";

  std::cout<<"   This program is distributed in the hope that it will be useful,        "<<"\n";
  std::cout<<"   but WITHOUT ANY WARRANTY; without even the implied warranty of         "<<"\n";
  std::cout<<"   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           "<<"\n";
  std::cout<<"   GNU General Public License for more details.                           "<<"\n";

  std::cout<<"   You should have received a copy of the GNU General Public License      "<<"\n";
  std::cout<<"   along with this program.  If not, see <http://www.gnu.org/licenses/>.  "<<"\n";

  std::cout<<"   For help with FERRET please contact J. Mangeri <john.mangeri@list.lu>  "<<"\n";
  std::cout<<"   and be sure to track new changes at github.com/mangerij/ferret         "<<"\n";
  std::cout<<"__________________________________________________________________________"<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"                                                                          "<<"\n";
  std::cout<<"  Initializing simulation:                                                "<<"\n";
  std::cout<<"__________________________________________________________________________"<<"\n";

  // Execute the application
  app->run();

  return 0;
}
