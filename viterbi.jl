# Author: Roman Sinayev
# License: BSD Style.
# Viterbi algorithm implementation and example written in Julia Language
#
# The idea is to decode the most probable order of a hidden markov model 
# based on indirect observations.

# In the example below, we attempt to determine where a person is based on
# their cell phone service strength.  We also have fake probabilities associated
# with each state(location) and with transitions between locations and with
# estimated signal strength probability at each location.

# In the given example, from the observations of 
# "three_bars", "five_bars","no_service","five_bars", "three_bars"
# we determine that the person goes from home to the street to starbucks, to
# the street again and then to work.  Such a sequence may not be immediately
# apparent from just observing the probabilities

# such information may be useful for apps that want to locate a person or try
# to determine if a person has visited a specific location for marketing
# purposes

states = ("home", "street", "starbucks", "work")
obs = ("three_bars", "five_bars","no_service","five_bars", "three_bars")
start_probs = ["home"=>0.9, "street"=>0.1, "starbucks"=>0.0, "work"=>0.0]
trans_probs = [ "home"=> 
                  ["home"=>0.4,"street"=>0.6,"starbucks"=>0.0,"work"=>0.0],
                "street"=>
                  ["home"=>0.0,"street"=>0.2,"starbucks"=>0.4,"work"=>0.4],
                "starbucks"=>
                  ["home"=>0.0,"street"=>0.5,"starbucks"=>0.5,"work"=>0.0],
                "work"=>
                  ["home"=>0.0,"street"=>0.1,"starbucks"=>0.0,"work"=>0.9]
              ]
emis_probs = ["home"=>
                ["five_bars"=>0.1, "no_service"=>0.3, "three_bars"=>0.6],
              "street"=>
                ["five_bars"=>0.8, "no_service"=>0.0, "three_bars"=>0.2],
              "starbucks"=>
                ["five_bars"=>0.1, "no_service"=>0.8, "three_bars"=>0.1],
              "work"=>
                ["five_bars"=>0.1, "no_service"=>0.0, "three_bars"=>0.9]
              ]

function viterbi(obs, states, start_probs, trans_probs, emis_probs)
  V = Array(Dict,0)
  # V is a series of dictionaries with location linked to probability of
  # that location at each time point
  push!(V, [state=>start_probs[state] * emis_probs[state][obs[1]] 
                                                    for state in states])
  path = Dict{ASCIIString,Array{ASCIIString,1}}()
  # path is the traveled path with dictionary pointing to the last item in
  # the list for easy backtracking
  for state in states
    path[state] = [state]
  end
  for t=2:length(obs)
    new_path = Dict{ASCIIString,Array{ASCIIString,1}}()
    for state in states
      prob, old_state = max(
        [(V[t-1][s] * trans_probs[s][state] * emis_probs[state][obs[t]], s) 
        for s in states]
        )
      # here we determined the best current location from the O(s^2) state
      # combinations
      push!(V, Dict{ASCIIString,Float64}())
      V[t][state] = prob
      new_path[state] = cat(1, path[old_state], [state])
    end
    path = new_path
  end
  prob, state = max([(V[length(obs)][s], s) for s in states])
  return prob, path[state]
end

prob,path = viterbi(obs, states, start_probs, trans_probs, emis_probs)

println(
  "the final estimated path is\n$(path)\nwith probability of $(round(prob,3))")