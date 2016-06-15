vcontext("logos")
data(sequences)

ggplot(data = ggfortify(sequences, "peptide")) +      
  geom_logo(aes(x=position, y=bits, group=element, 
     label=element, fill=interaction(Polarity, Water)),
     alpha = 0.6)  +
  scale_fill_brewer(palette="Paired") +
  theme(legend.position = "bottom")
save_vtest(desc="first plot")

ggplot(data = ggfortify(sequences, "peptide", treatment = "class")) + 
  geom_logo(aes(x=class, y=bits, group=element, 
     label=element, fill=element)) + 
  facet_wrap(~position, ncol=18) +
  theme(legend.position = "bottom")
save_vtest(desc="facet by position")

ggplot(data = ggfortify(sequences, "peptide", treatment = "class")) + 
  geom_logo(aes(x=position, y=bits, group=element, label=element, fill=element)) + 
  facet_wrap(~class, ncol=1) + theme_bw()
save_vtest(desc="facet by class")


ggplot(data = ggfortify(sequences, "peptide", treatment = "class")) + 
  geom_logo(aes(x=class, y=bits, group=element, 
                label=element, fill=interaction(Polarity, Water))) + 
  scale_fill_brewer("Amino-acids properties", palette="Paired") + 
  facet_wrap(~position, ncol=18) + 
  theme(legend.position="bottom") + 
  xlab("") + ylab("Shannon information in bits")
save_vtest(desc="facet by position, adjust color & label")

end_vcontext()