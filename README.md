# How to use

Ebi is a command line tool that requires neither installation nor internet access.

# Getting started with development

1. Install Rustup
    https://www.rust-lang.org/tools/install

1. Log out and in again

1. Install Visual Studio Code

1. Install extension 'rust-analyzer' in Visual Studio Code
    - https://code.visualstudio.com/docs/languages/rust

1. To run Ebi, use the terminal of Visual Studio Code to give the command "cargo run --" instead of "Ebi". Everything else is equivalent to the commands mentioned in the manual.

1. The SLPN model used for case study is called permit2018_sm.slpn. Given trace "Permit submitted by employee, Permit final approved by supervisor, Start trip, Declaration submitted by employee, Declaration final approved by supervisor, Request payment, End trip, Payment handled", the following command can be used to compute a stochastic alignment with balance factor set to 1:
```
   cargo run probability exptra testfiles/permit2018_sm.slpn 1 "Permit SUBMITTED by EMPLOYEE" "Permit FINAL_APPROVED by SUPERVISOR" "Start trip" "Declaration SUBMITTED by EMPLOYEE" "Declaration FINAL_APPROVED by SUPERVISOR" "Request Payment" "End trip" "Payment Handled"
```
The following command can be used to compute a stochastic alignment with balance factor set to 0.5:
```
   cargo run probability exptra testfiles/permit2018_sm.slpn 0.5 "Permit SUBMITTED by EMPLOYEE" "Permit FINAL_APPROVED by SUPERVISOR" "Start trip" "Declaration SUBMITTED by EMPLOYEE" "Declaration FINAL_APPROVED by SUPERVISOR" "Request Payment" "End trip" "Payment Handled"
```
The probability of the identified model paths m1 = <Permit submitted by employee, Permit approved by administration, Permit final approved by supervisor, Start trip, Declaration submitted by employee, Declaration approved by administration, Declaration final approved by supervisor, Request payment, Payment handled> can be queried using the following cmd:

```
cargo run probability trac testfiles/permit2018_sm.slpn "Permit SUBMITTED by EMPLOYEE" "Permit APPROVED by ADMINISTRATION" "Permit FINAL_APPROVED by SUPERVISOR" "Start trip" "Declaration SUBMITTED by EMPLOYEE" "Declaration APPROVED by ADMINISTRATION" "Declaration FINAL_APPROVED by SUPERVISOR" "Request Payment" "Payment Handled"
```

Likewise, the probability of the identified model paths m0.5= <Permit submitted by employee, Permit approved by administration, Permit final approved by supervisor, Start trip, Declaration submitted by employee, Declaration approved by administration, Declaration final approved by supervisor, Request payment, Payment handled> can be queried using the following cmd:
```
cargo run probability trace testfiles/permit2018_sm.slpn "Permit SUBMITTED by EMPLOYEE" "Permit APPROVED by ADMINISTRATION" "Permit FINAL_APPROVED by SUPERVISOR" "Start trip" "End trip" "Declaration SUBMITTED by EMPLOYEE" "Declaration APPROVED by ADMINISTRATION" "Declaration FINAL_APPROVED by SUPERVISOR" "Request Payment" "Payment Handled"
```
