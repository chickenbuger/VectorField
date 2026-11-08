#ifndef FIELDLAB_FRAME_SCHEDULER_H
#define FIELDLAB_FRAME_SCHEDULER_H
#include <algorithm>
#include <cmath>
namespace fieldlab {
// Interactive work is bounded; captures keep the complete fixed-step sequence.
// A physics step cannot be interrupted: one slow step can exceed the CPU budget.
class FrameScheduler {
public:
    static double FixedStep() { return 1.0/120.0; }
    void Reset() { pending_=discarded_=0; }
    void AddElapsed(double elapsed,bool interactive) {
        if(!(elapsed>0)||!std::isfinite(elapsed))return;
        const double accepted=(std::min)(0.1,elapsed);
        discarded_+=elapsed-accepted;pending_+=accepted;
        if(interactive){
            const double limit=4*FixedStep();
            if(pending_>limit){discarded_+=pending_-limit;pending_=limit;}
        }
    }
    bool ShouldStep(int completed,double cpuMs,bool interactive) const {
        return pending_>=FixedStep() && (!interactive ||
            (completed<3 && (completed==0 || cpuMs<12.0)));
    }
    void ConsumeStep() { pending_=(std::max)(0.0,pending_-FixedStep()); }
    double PendingSeconds() const { return pending_; }
    double DiscardedSeconds() const { return discarded_; }
private:
    double pending_=0,discarded_=0;
};
}
#endif
