///////////////////////////////////////////////////////////////////////////
// Berenger Bramas INRIA - 2014
// Code provided under GNU Lesser General Public License
//
//
// This file ask the cpuid to get access to CPU properties.
// The file contains 3 mains parts:
// × First part is a wrapper in case we are on Windows or Linux
// × Second part is several call to the function and fill of an list
// × Third is out of scope, it prints the state and the properties
//   in a strict format in order to post process
///////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// Part 1:
// Defines cpuid:
// × A Wrapper if we are on windows
// × A call to assembly else
///////////////////////////////////////////////////////////////////////////

enum RegistersNum {
    EaxRegister = 0,
    EbxRegister,
    EcxRegister,
    EdxRegister
};

#ifdef _WIN32

// On windows __cpuid exists: http://msdn.microsoft.com/en-us/library/hskdteyh(v=vs.90).aspx
// void __cpuid(int CPUInfo[4],int InfoType);
// we would like to have the same name for not windows
#define cpuid    __cpuid

#else

// Else we have to ask the CPU directly by executin cpuid.
// eax should contains the information querry argument.
// Then we have to take the results from the different registers.
//
//    From : http://www.ibiblio.org/gferg/ldp/GCC-Inline-Assembly-HOWTO.html
//
//    asm ( assembler template
//        : output operands                  // optional
//        : input operands                   // optional
//        : list of clobbered registers      // optional
//        );
//
//    +---+--------------------+
//    | r |    Register(s)     |
//    +---+--------------------+
//    | a |   %eax, %ax, %al   |
//    | b |   %ebx, %bx, %bl   |
//    | c |   %ecx, %cx, %cl   |
//    | d |   %edx, %dx, %dl   |
//    | S |   %esi, %si        |
//    | D |   %edi, %di        |
//    +---+--------------------+
//


//  GCC Inline Assembly but with the same prototype as windows
void cpuid(int CPUInfo[4],int InfoType){
    __asm__ __volatile__ (
        "cpuid":            // Execute this instruction
        "=a" (CPUInfo[EaxRegister]),  // Store eax in 0
        "=b" (CPUInfo[EbxRegister]),  // Store ebx in 1
        "=c" (CPUInfo[EcxRegister]),  // Store ecx in 2
        "=d" (CPUInfo[EdxRegister]) : // Store edx in 3
        "a" (InfoType)      // Input InfoType in eax before instruction
    );
}

#endif


bool CPUInfoGetEAX(const int CPUInfo[4], const int position){
    return (CPUInfo[EaxRegister] & ((int)1 << position)) != 0;
}

bool CPUInfoGetEBX(const int CPUInfo[4], const int position){
    return (CPUInfo[EbxRegister] & ((int)1 << position)) != 0;
}

bool CPUInfoGetECX(const int CPUInfo[4], const int position){
    return (CPUInfo[EcxRegister] & ((int)1 << position)) != 0;
}

bool CPUInfoGetEDX(const int CPUInfo[4], const int position){
    return (CPUInfo[EdxRegister] & ((int)1 << position)) != 0;
}

///////////////////////////////////////////////////////////////////////////
// Part 2:
// Call the cpuid function and ask for particular information.
// In our case we want to use these information to print it (and later use it
// in a CMake file).
// So you can change this file to get more informations and do something else with them.
//
//    From 64-ia-32-architectures-software-developer-instruction-set-reference-manual-325383.pdf
//    Or recommanded : AMD CPUID_Specification.pdf
//    We know in part CPUID—CPU Identification that a call to the cpuid instruction fill the registers
//    with the cpu property.
///////////////////////////////////////////////////////////////////////////

#include <string>
#include <list>

struct CpuProperty {
    CpuProperty(const char inName[], const bool IsEnable)
        : name(inName), enabled(IsEnable){
    }

    std::string name;
    bool enabled;
};

std::list<CpuProperty> getProperties(){
    std::list<CpuProperty> properties;

    // To store the registers value
	int info[4];

    // Basic CPUID Information
    cpuid(info, 0);
    // The largest CPUID standard-function input value supported by the processor implementation.
    const int limitStandardFunction = info[EaxRegister];

    // Extended Function CPUID Information
    cpuid(info, 0x80000000);
    // The largest CPUID extended-function input value supported by the processor implementation
    int limitExtendedFunction = info[EaxRegister];

	//  Detect Instruction Set
    if (limitStandardFunction >= 1){
        cpuid(info,0x00000001); // Basic CPUID Information
        /*
        0x00000001 - EDX :
            31:29 Reserved.
            28 HTT: Hyper-Threading Technology. Indicates either that there is more than one thread per CPU core
              or more than one CPU core per processor. AMD currently does not support more than one thread per
             CPU core. See “Legacy Method” on page 23.
            27 Reserved.
            26 SSE2: SSE2 extensions. See Appendix D “CPUID Feature Sets” in APM3.
            25 SSE: SSE extensions. See Appendix D “CPUID Feature Sets” in APM3 appendix and “64-Bit Media
                    Programming” in APM1.
            24 FXSR: FXSAVE and FXRSTOR instructions. See “FXSAVE” and “FXRSTOR” in APM4.
            23 MMX: MMXTM instructions. See Appendix D “CPUID Feature Sets” in APM3 and “128-Bit Media
                      and Scientific Programming” in APM1.
            22:20 Reserved.
            19 18 Reserved.
            17 PSE36: Page-size extensions. The PDE[20:13] supplies physical address [39:32]. See “Page Translation and Protection” in APM2.
            16 PAT: Page attribute table. PCD, PWT, and PATi are used to alter memory type. See “Page-Attribute Table Mechanism” in APM2.
            15 CMOV: Conditional move instructions, CMOV, FCMOV. See “CMOV”, “FCMOV” in APM3.
            14 MCA: Machine check architecture, MCG_CAP. See “Machine Check Mechanism” in APM2.
            13 PGE: Page global extension, CR4.PGE. See “Page Translation and Protection” in APM2.
            12 MTRR: Memory-type range registers. MTRRcap supported. See “Page Translation and Protection” in APM2.
            11 SysEnterSysExit: SYSENTER and SYSEXIT instructions. See “SYSENTER”, “SYSEXIT“ in APM3.
            10 Reserved
            9 APIC. Advanced programmable interrupt controller (APIC) exists and is enabled. See “Exceptions
               and Interrupts” in APM2.
            8 CMPXCHG8B: CMPXCHG8B instruction. See “CMPXCHG8B” in APM3.
            7 MCE: Machine check exception, CR4.MCE. See “Machine Check Mechanism” in APM2.
            6 PAE: Physical-address extensions (PAE), support for physical addresses ≥ 32b. Number of physical
               address bits above 32b is implementation specific. See “Page Translation and Protection” in APM2.
            5 MSR: AMD model-specific registers (MSRs), with RDMSR and WRMSR instructions. See “Model
               Specific Registers” in APM2.
            4 TSC: Time stamp counter. RDTSC and RDTSCP instruction support. See “Debug and Performance
               Resources” in APM2.
            3 PSE: Page-size extensions (4 MB pages). See “Page Translation and Protection” in APM2.
            2 DE: Debugging extensions, I/O breakpoints, CR4.DE. See “Debug and Performance Resources” in
                 APM2.
            1 VME: Virtual-mode enhancements, CR4.VME, CR4.PVI, software interrupt indirection, expansion
             of the TSS with the software, indirection bitmap, EFLAGS.VIF, EFLAGS.VIP. See “System
              Resources” in APM2.
            0 FPU: x87 floating point unit on-chip. See “x87 Floating Point Programming” in APM1
         */

        properties.push_back(CpuProperty("MMX", CPUInfoGetEDX(info, 23)));
        properties.push_back(CpuProperty("SSE", CPUInfoGetEDX(info, 25)));
        properties.push_back(CpuProperty("SSE2", CPUInfoGetEDX(info, 26)));

        /*
         0x00000001 - ECX :
            0 SSE3 Streaming SIMD Extensions 3 (SSE3). A value of 1 indicates the processor supports this
            technology.
            1 PCLMULQDQ PCLMULQDQ. A value of 1 indicates the processor supports the PCLMULQDQ instruction
            2 DTES64 64-bit DS Area. A value of 1 indicates the processor supports DS area using 64-bit layout
            3 MONITOR MONITOR/MWAIT. A value of 1 indicates the processor supports this feature.
            4 DS-CPL CPL Qualified Debug Store. A value of 1 indicates the processor supports the extensions to the
            Debug Store feature to allow for branch message storage qualified by CPL.
            5 VMX Virtual Machine Extensions. A value of 1 indicates that the processor supports this technology
            6 SMX Safer Mode Extensions. A value of 1 indicates that the processor supports this technology. See
            Chapter 5, “Safer Mode Extensions Reference”.
            7 EIST Enhanced Intel SpeedStep® technology. A value of 1 indicates that the processor supports this
            technology.
            8 TM2 Thermal Monitor 2. A value of 1 indicates whether the processor supports this technology.
            9 SSSE3 A value of 1 indicates the presence of the Supplemental Streaming SIMD Extensions 3 (SSSE3). A
            value of 0 indicates the instruction extensions are not present in the processor
            10 CNXT-ID L1 Context ID. A value of 1 indicates the L1 data cache mode can be set to either adaptive mode
            or shared mode. A value of 0 indicates this feature is not supported. See definition of the
            IA32_MISC_ENABLE MSR Bit 24 (L1 Data Cache Context Mode) for details.
            11 SDBG A value of 1 indicates the processor supports IA32_DEBUG_INTERFACE MSR for silicon debug.
            12 FMA A value of 1 indicates the processor supports FMA extensions using YMM state.
            13 CMPXCHG16B CMPXCHG16B Available. A value of 1 indicates that the feature is available. See the
            “CMPXCHG8B/CMPXCHG16B—Compare and Exchange Bytes” section in this chapter for a
            description.
            14 xTPR Update
            Control
            xTPR Update Control. A value of 1 indicates that the processor supports changing
            IA32_MISC_ENABLE[bit 23].
            15 PDCM Perfmon and Debug Capability: A value of 1 indicates the processor supports the performance
            and debug feature indication MSR IA32_PERF_CAPABILITIES.
            16 Reserved Reserved
            17 PCID Process-context identifiers. A value of 1 indicates that the processor supports PCIDs and that
            software may set CR4.PCIDE to 1.
            18 DCA A value of 1 indicates the processor supports the ability to prefetch data from a memory mapped
            device.
            19 SSE4.1 A value of 1 indicates that the processor supports SSE4.1.
            20 SSE4.2 A value of 1 indicates that the processor supports SSE4.2.
            21 x2APIC A value of 1 indicates that the processor supports x2APIC feature.
            22 MOVBE A value of 1 indicates that the processor supports MOVBE instruction.
            23 POPCNT A value of 1 indicates that the processor supports the POPCNT instruction.
            24 TSC-Deadline A value of 1 indicates that the processor’s local APIC timer supports one-shot operation using a
            TSC deadline value.
            25 AESNI A value of 1 indicates that the processor supports the AESNI instruction extensions.
            26 XSAVE A value of 1 indicates that the processor supports the XSAVE/XRSTOR processor extended states
            feature, the XSETBV/XGETBV instructions, and XCR0.
            27 OSXSAVE A value of 1 indicates that the OS has set CR4.OSXSAVE[bit 18] to enable the XSAVE feature set.
            28 AVX A value of 1 indicates the processor supports the AVX instruction extensions.
            29 F16C A value of 1 indicates that processor supports 16-bit floating-point conversion instructions.
            30 RDRAND A value of 1 indicates that processor supports RDRAND instruction.
            31 Not Used Always returns 0.
          */
        properties.push_back(CpuProperty("SSE3", CPUInfoGetECX(info,  0)));
        properties.push_back(CpuProperty("SSSE3", CPUInfoGetECX(info,  9)));
        properties.push_back(CpuProperty("SSE41", CPUInfoGetECX(info, 19)));
        properties.push_back(CpuProperty("SSE42", CPUInfoGetECX(info, 20)));
        properties.push_back(CpuProperty("AVX",  CPUInfoGetECX(info, 28)));
        properties.push_back(CpuProperty("FMA3", CPUInfoGetECX(info, 12)));
	}

    if (limitExtendedFunction >= 0x80000001){
        cpuid(info,0x80000001); // Extended Function CPUID Information
        /*
        0x80000001 - EDX :
            31 3DNow: 3DNow!TM instructions. See Appendix D “Instruction Subsets and CPUID Feature Sets” in APM3.
            30 3DNowExt: AMD extensions to 3DNow! instructions. See Appendix D “Instruction Subsets and
                CPUID Feature Sets” in APM3.
            29 LM: Long mode. See “Processor Initialization and Long-Mode Activation” in APM2.
            28 Reserved.
            27 RDTSCP: RDTSCP instruction. See “RDTSCP” in APM3.
            26 Page1GB: 1-GB large page support. See “1-GB Paging Support” in APM2.
            25 FFXSR: FXSAVE and FXRSTOR instruction optimizations. See “FXSAVE” and “FXRSTOR” in APM4.
            24 FXSR: FXSAVE and FXRSTOR instructions. Same as CPUID Fn0000_0001_EDX[FXSR].
            23 MMX: MMXTM instructions. Same as CPUID Fn0000_0001_EDX[MMX].
            22 MmxExt: AMD extensions to MMX instructions. See Appendix D “Instruction Subsets and CPUID
                Feature Sets” in APM3 and “128-Bit Media and Scientific Programming” in APM1.
            21 Reserved.
            20 NX: No-execute page protection. See “Page Translation and Protection” in APM2.
            19:18 Reserved.
            17 PSE36: Page-size extensions. Same as CPUID Fn0000_0001_EDX[PSE36].
            16 PAT: Page attribute table. Same as CPUID Fn0000_0001_EDX[PAT].
            15 CMOV: Conditional move instructions. Same as CPUID Fn0000_0001_EDX[CMOV]
            14 MCA: Machine check architecture. Same as CPUID Fn0000_0001_EDX[MCA].
            13 PGE: Page global extension. Same as CPUID Fn0000_0001_EDX[PGE].
            12 MTRR: Memory-type range registers. Same as CPUID Fn0000_0001_EDX[MTRR].
            11 SysCallSysRet: SYSCALL and SYSRET instructions. See “SYSCALL” and “SYSRET” in APM3.
            10 Reserved.
            9 APIC. Advanced programmable interrupt controller. Same as CPUID Fn0000_0001_EDX[APIC].
            8 CMPXCHG8B: CMPXCHG8B instruction. Same as CPUID Fn0000_0001_EDX[CMPXCHG8B].
            7 MCE: Machine check exception. Same as CPUID Fn0000_0001_EDX[MCE].
            6 PAE: Physical-address extensions. Same as CPUID Fn0000_0001_EDX[PAE].
            5 MSR: AMD model-specific registers. Same as CPUID Fn0000_0001_EDX[MSR].
            4 TSC: Time stamp counter. Same as CPUID Fn0000_0001_EDX[TSC].
            3 PSE: Page-size extensions. Same as CPUID Fn0000_0001_EDX[PSE].
            2 DE: Debugging extensions. Same as CPUID Fn0000_0001_EDX[DE].
            1 VME: Virtual-mode enhancements. Same as CPUID Fn0000_0001_EDX[VME].
            0 FPU: x87 floating-point unit on-chip. Same as CPUID Fn0000_0001_EDX[FPU].
          */
        properties.push_back(CpuProperty("x64", CPUInfoGetEDX(info, 29)));
        /*
        0x80000001 - ECX :
            31:14 Reserved.
            13 WDT: Watchdog timer support.
            12 SKINIT: SKINIT, STGI, and DEV support.
            11:10 Reserved.
            9 OSVW: OS visible workaround. Indicates OS-visible workaround support. See “OS Visible Work-
             around (OSVW) Information” in APM2.
            8 3DNowPrefetch: PREFETCH and PREFETCHW instruction support. See “PREFETCH” and “PREFETCHW” in APM3.
            7 MisAlignSse: Misaligned SSE mode. See “Misaligned Access Support Added for SSE Instructions”
                 in APM1.
            6 SSE4A: EXTRQ, INSERTQ, MOVNTSS, and MOVNTSD instruction support. See “EXTRQ”,
                 “INSERTQ”, “MOVNTSS”, and “MOVNTSD” in APM4.
            5 ABM: Advanced bit manipulation. LZCNT instruction support. See “LZCNT” in APM3.
            4 AltMovCr8: LOCK MOV CR0 means MOV CR8. See “MOV(CRn)” in APM3.
            3 ExtApicSpace: This bit indicates the presence of extended APIC register space starting at offset
             400h from the “APIC Base Address Register,” as specified in the BKDG..
            2 SVM: Secure virtual machine feature. See “Secure Virtual Machine” in APM2.
            1 CmpLegacy: Core multi-processing legacy mode. See “Legacy Method” on page 23.
            0 LahfSahf: LAHF and SAHF instruction support in 64-bit mode. See “LAHF” and “SAHF” in   APM3.
          */
        properties.push_back(CpuProperty("SSE4a", CPUInfoGetECX(info,  6)));
        properties.push_back(CpuProperty("FMA4", CPUInfoGetECX(info, 16)));
        properties.push_back(CpuProperty("XOP", CPUInfoGetECX(info, 11)));
	}

    return properties;
}


///////////////////////////////////////////////////////////////////////////
// Part 3:
// Print the information in a format to use it with CMake
///////////////////////////////////////////////////////////////////////////

#include <iostream>

int main(){
    const std::list<CpuProperty> properties = getProperties();

    const std::list<CpuProperty>::const_iterator endIterProperties = properties.end();
    for(std::list<CpuProperty>::const_iterator iterProperties = properties.begin()
            ; iterProperties != endIterProperties
            ; ++iterProperties){
        // Print the status
        std::cout << (*iterProperties).name << "=" << ((*iterProperties).enabled?"TRUE":"FALSE") << ";";
    }

    return 0;
}
