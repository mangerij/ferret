# ferret
Open-source MOOSE application "Ferret" for parallel mesoscale simulations of ferroic and related electronic materials

See the website at https://mangerij.github.io/ferret/ for more information.

# note on LLM-inspired code changes

Proposed code changes to the Ferret repository using LLM tools are welcome.
However, it should be stressed that existing Ferret objects should be inspected that may already cover your use case before new ones are invented.
These tools often duplicate existing functionality that may provide subtle solution differences that do not pass benchmark scrutiny.

New code should include tests that pass the CI, which runs against the underlying MOOSE framework. In general, this is a good practice for research-grade use of Ferret.

As of 23-09-26, legacy objects have been carefully marked (and some removed) which may have (previous to-this-date) caused confusion for LLMs reading this repository.
Tests also have been updated to cover the synthesized code branch that covers a variety of examples for cubic perovskites (e.g., thin films under epitaxy).

More changes are coming with regards to magnetic (and multiferroic) systems.

If you are unsure of your changes, please contact us. We are happy to help.
