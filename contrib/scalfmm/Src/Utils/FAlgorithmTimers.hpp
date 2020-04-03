// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

#ifndef FALGORITHMTIMERS_HPP
#define FALGORITHMTIMERS_HPP

/**
 * @brief Collection of timers for FMM operators.
 *
 * This class provide a way for the different algorithms to
 * store the time spent in each operator.
 */
class FAlgorithmTimers{
public:
    /// The timer names
    enum FTimers {
        P2MTimer,
        M2MTimer,
        M2LTimer,
        L2LTimer,
        L2PTimer,
        P2PTimer,
        NearTimer,
        nbTimers   ///< Timer count
    };

protected:
    /// Timer array
    FTic Timers[nbTimers];

public:
    /// Constructor: resets all timers
    FAlgorithmTimers()
    {
        for(int i = 0; i < nbTimers ; ++i){
            Timers[i].reset();
        }
    }

    /// Default copy contructor
    FAlgorithmTimers(const FAlgorithmTimers&) = default;
    /// Default move contructor
    FAlgorithmTimers(FAlgorithmTimers&&) = default;

    /// Returns the timer array
    const FTic * getAllTimers() const {
        return Timers;
    }

    /// Returns the timer count
    int getNbOfTimerRecorded() const {
        return nbTimers;
    }

    /// Elapsed time between last FTic::tic() and FTic::tac() for given timer.
    double getTime(FTimers OpeTimer) const{
        //assert to verify size
        return Timers[OpeTimer].elapsed();
    }

    /// Cumulated time between all FTic::tic() and FTic::tac() for given timer.
    double getCumulatedTime(FTimers OpeTimer) const{
        //assert to verify size
        return Timers[OpeTimer].cumulated();
    }

};

#endif // FALGORITHMTIMERS_HPP
