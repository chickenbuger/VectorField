#ifndef FIELDLAB_FRAME_SCHEDULER_TESTS_H
#define FIELDLAB_FRAME_SCHEDULER_TESTS_H
#include "../opgnel_project_smoke/FrameScheduler.h"
#include <iostream>
#include <limits>
namespace fieldlab {
inline int RunFrameSchedulerTests() {
    int passed=0,failed=0;
    const auto check=[&](const char* name,bool ok){(ok?passed:failed)++;std::cout<<(ok?"PASS ":"FAIL ")<<name<<'\n';};
    FrameScheduler s;s.AddElapsed(0.1,true);
    check("interactive backlog is bounded with discarded wall time accounted",std::abs(s.PendingSeconds()-4.0/120)<1e-12&&std::abs(s.DiscardedSeconds()-(0.1-4.0/120))<1e-12);
    check("interactive callback permits its first complete step",s.ShouldStep(0,100,true));s.ConsumeStep();
    check("interactive CPU budget returns control before another expensive step",!s.ShouldStep(1,12,true)&&s.ShouldStep(1,11.9,true));
    check("interactive callback has a three-step upper bound",!s.ShouldStep(3,0,true)&&s.ShouldStep(2,0,true));
    s.Reset();s.AddElapsed(0.1,false);int n=0;while(s.ShouldStep(n,100,false)){s.ConsumeStep();++n;}
    check("capture scheduling retains all twelve fixed steps",n==12&&s.DiscardedSeconds()==0);
    s.Reset();s.AddElapsed(1.0/240,true);check("fractional frame time waits for a full physical step",!s.ShouldStep(0,0,true));s.AddElapsed(1.0/240,true);
    check("fractional frame times combine without variable dt",s.ShouldStep(0,0,true));s.ConsumeStep();
    check("consuming a fixed step preserves zero remainder",s.PendingSeconds()<1e-12);
    s.AddElapsed(-1,true);s.AddElapsed(std::numeric_limits<double>::infinity(),true);
    check("invalid wall-clock intervals cannot create simulation debt",s.PendingSeconds()<1e-12);
    s.AddElapsed(0.5,false);check("wall-clock stall clamp reports discarded time",std::abs(s.DiscardedSeconds()-0.4)<1e-12);
    s.Reset();check("scene reset clears pending and discarded wall time",s.PendingSeconds()==0&&s.DiscardedSeconds()==0);
    std::cout<<"FrameScheduler tests: "<<passed<<" passed, "<<failed<<" failed, "<<passed+failed<<" total.\n";
    return failed?1:0;
}
}
#endif
